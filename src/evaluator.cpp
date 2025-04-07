#include "evaluator.h"
#include "structures.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream> // For placeholder messages
#include <set>      // Needed for std::set

// --- Constructor ---

Evaluator::Evaluator(Circuit &circuit, const std::map<std::string, Point> &placement)
    : circuit_ref(circuit), placement_ref(placement)
{
    // Initialization if needed
}

// --- Core Evaluation Metrics ---

double Evaluator::calculate_total_hpwl()
{
    double total_hpwl = 0.0;
    for (const auto &net_pair : circuit_ref.net_map)
    {
        total_hpwl += calculate_net_hpwl(net_pair.second);
    }
    return total_hpwl;
}

double Evaluator::calculate_placement_area()
{
    if (placement_ref.empty())
    {
        return 0.0;
    }

    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    for (const auto &pair : placement_ref)
    {
        // Only consider movable nodes for the placement area calculation
        auto node_it = circuit_ref.cell_map.find(pair.first);
        if (node_it == circuit_ref.cell_map.end())
        {
            std::cerr << "Warning: Node ID '" << pair.first << "' found in placement but not in circuit cell_map." << std::endl;
            continue; // Skip if node definition not found
        }
        const Node *node = &(node_it->second); // Get pointer to the node in the map

        if (node && node->type == MOVABLE)
        {
            const Point &pos = pair.second;
            // Access width and height directly from Node struct
            const double node_width = node->width;
            const double node_height = node->height;

            min_x = std::min(min_x, pos.x);
            max_x = std::max(max_x, pos.x + node_width); // Use node width
            min_y = std::min(min_y, pos.y);
            max_y = std::max(max_y, pos.y + node_height); // Use node height
        }
    }

    // If no movable nodes were placed, area is 0
    if (min_x == std::numeric_limits<double>::max())
    {
        return 0.0;
    }

    return (max_x - min_x) * (max_y - min_y);
}

double Evaluator::estimate_total_overflow()
{
    if (!demand_estimated)
    {
        estimate_routing_demand();
    }

    double total_overflow = 0.0;
    // Calculate overflow for each routing grid edge
    for (const auto &edge : routing_demand)
    {
        if (edge.second > edge.first.second) // If demand exceeds capacity
        {
            total_overflow += (edge.second - edge.first.second);
        }
    }
    return total_overflow;
}

double Evaluator::estimate_critical_delay()
{
    double max_delay = 0.0;

    // Calculate delay for each net
    for (const auto &net_pair : circuit_ref.net_map)
    {
        const Net &net = net_pair.second;
        double net_delay = calculate_elmore_delay(net);
        max_delay = std::max(max_delay, net_delay);
    }

    return max_delay;
}

// --- Feedback Generation ---

std::vector<std::string> Evaluator::get_congested_regions_feedback()
{
    // Placeholder implementation
    std::cout << "Warning: get_congested_regions_feedback() is not fully implemented." << std::endl;
    if (!demand_estimated)
    {
        estimate_routing_demand();
    }
    // A real implementation would analyze the estimated routing demand vs. capacity
    // on the grid and identify regions exceeding a threshold.
    return {}; // Return empty vector for now
}

// --- Private Helper Methods ---

double Evaluator::calculate_net_hpwl(const Net &net)
{
    // Use net.nodes which contains Node* instead of non-existent net.pins
    if (net.nodes.empty())
    {
        return 0.0;
    }

    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    bool found_valid_center = false;

    // Iterate through the nodes connected by the net
    for (const Node *node_ptr : net.nodes)
    {
        if (!node_ptr)
        {
            std::cerr << "Warning: Null node pointer encountered in net '" << net.name << "' for HPWL calculation." << std::endl; // Use net.name
            continue;
        }

        const std::string &node_id = node_ptr->name; // Get node name

        // Node definition should exist since we got the pointer from Circuit structure,
        // but placement might be missing for movable nodes.
        const Node *node = node_ptr; // Use the pointer directly

        Point node_pos;
        // Get position from placement map if movable, or use fixed position from node->pos
        if (node->type == MOVABLE) // Use node->type
        {
            auto it = placement_ref.find(node_id);
            if (it == placement_ref.end())
            {
                // This should ideally not happen if placement is complete for movable nodes
                std::cerr << "Warning: Movable node '" << node_id << "' not found in placement map for net '" << net.name << "' HPWL calculation." << std::endl; // Use net.name
                continue;                                                                                                                                        // Skip this node if it isn't placed
            }
            node_pos = it->second;
        }
        else // Handles FIXED, TERMINAL, TERMINAL_NI
        {
            node_pos = node->pos; // Use the fixed position from the Node struct (node->pos)
        }

        // Calculate the center of the node using node->width and node->height
        double center_x = node_pos.x + node->width / 2.0;
        double center_y = node_pos.y + node->height / 2.0;

        min_x = std::min(min_x, center_x);
        max_x = std::max(max_x, center_x);
        min_y = std::min(min_y, center_y);
        max_y = std::max(max_y, center_y);
        found_valid_center = true;
    }

    // If no valid node centers were found (e.g., all nodes missing from placement or circuit)
    if (!found_valid_center)
    {
        return 0.0;
    }

    return (max_x - min_x) + (max_y - min_y);
}

void Evaluator::estimate_routing_demand()
{
    // Clear previous demand
    routing_demand.clear();

    // Define routing grid parameters
    const double GRID_SIZE = 1.0;    // Size of each routing grid cell
    const double WIRE_SPACING = 0.1; // Minimum spacing between wires

    // Calculate routing demand for each net
    for (const auto &net_pair : circuit_ref.net_map)
    {
        const Net &net = net_pair.second;
        if (net.nodes.empty())
            continue;

        // Calculate bounding box of the net
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::lowest();

        for (const Node *node : net.nodes)
        {
            if (!node)
                continue;

            Point pos = node->type == MOVABLE ? placement_ref.at(node->name) : node->pos;
            min_x = std::min(min_x, pos.x);
            max_x = std::max(max_x, pos.x + node->width);
            min_y = std::min(min_y, pos.y);
            max_y = std::max(max_y, pos.y + node->height);
        }

        // Calculate number of grid cells needed
        int grid_cells_x = static_cast<int>((max_x - min_x) / GRID_SIZE) + 1;
        int grid_cells_y = static_cast<int>((max_y - min_y) / GRID_SIZE) + 1;

        // Add demand to each grid cell in the bounding box
        for (int x = 0; x < grid_cells_x; ++x)
        {
            for (int y = 0; y < grid_cells_y; ++y)
            {
                std::pair<int, int> grid_pos(x, y);
                routing_demand[grid_pos] += 1.0 / (WIRE_SPACING * WIRE_SPACING);
            }
        }
    }

    demand_estimated = true;
}

double Evaluator::calculate_elmore_delay(const Net &net)
{
    if (net.nodes.empty())
        return 0.0;

    // Find driver node (assuming first node is driver)
    const Node *driver = net.nodes[0];
    if (!driver)
        return 0.0;

    double total_delay = 0.0;

    // Calculate RC delay for each sink
    for (size_t i = 1; i < net.nodes.size(); ++i)
    {
        const Node *sink = net.nodes[i];
        if (!sink)
            continue;

        // Get positions from placement
        Point driver_pos = driver->type == MOVABLE ? placement_ref.at(driver->name) : driver->pos;
        Point sink_pos = sink->type == MOVABLE ? placement_ref.at(sink->name) : sink->pos;

        // Calculate Manhattan distance for wire length
        double wire_length = std::abs(driver_pos.x - sink_pos.x) + std::abs(driver_pos.y - sink_pos.y);

        // Simple RC delay model
        // R = resistance per unit length * length
        // C = capacitance per unit length * length
        const double R_PER_UNIT = 0.1; // Resistance per unit length
        const double C_PER_UNIT = 0.2; // Capacitance per unit length
        double R = R_PER_UNIT * wire_length;
        double C = C_PER_UNIT * wire_length;

        // Add driver resistance and sink capacitance
        // R += driver->resistance;
        // C += sink->capacitance;

        // Elmore delay = R * C
        total_delay = std::max(total_delay, R * C);
    }

    return total_delay;
}