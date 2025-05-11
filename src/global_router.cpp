#include "../include/global_router.h" // Adjust path as necessary
#include <cmath>                      // For std::ceil
#include <stdexcept>                  // For std::out_of_range

// RoutingGrid Constructor
InternalRoutingGrid::InternalRoutingGrid(int num_x_tiles, int num_y_tiles, int num_layers)
    : num_x_tiles_(num_x_tiles), num_y_tiles_(num_y_tiles), num_layers_(num_layers)
{
    if (num_x_tiles <= 0 || num_y_tiles <= 0 || num_layers <= 0)
    {
        // Or handle this error more gracefully, maybe throw an exception
        // For now, just ensuring grid_cells_ doesn't get invalid dimensions.
        // This case should ideally be caught by build_routing_grid before construction.
        this->num_x_tiles_ = 0;
        this->num_y_tiles_ = 0;
        this->num_layers_ = 0;
        return;
    }

    grid_cells_.resize(num_layers_);
    for (int l = 0; l < num_layers_; ++l)
    {
        grid_cells_[l].resize(num_x_tiles_);
        for (int i = 0; i < num_x_tiles_; ++i)
        {
            grid_cells_[l][i].resize(num_y_tiles_);
        }
    }
}

// Accessor for modifiable GridEdge data
GridEdge &InternalRoutingGrid::get_grid_cell_data(int layer, int x, int y)
{
    if (layer < 0 || layer >= num_layers_ ||
        x < 0 || x >= num_x_tiles_ ||
        y < 0 || y >= num_y_tiles_)
    {
        throw std::out_of_range("Grid coordinates out of range");
    }
    return grid_cells_[layer][x][y];
}

// Accessor for const GridEdge data
const GridEdge &InternalRoutingGrid::get_grid_cell_data(int layer, int x, int y) const
{
    if (layer < 0 || layer >= num_layers_ ||
        x < 0 || x >= num_x_tiles_ ||
        y < 0 || y >= num_y_tiles_)
    {
        throw std::out_of_range("Grid coordinates out of range");
    }
    return grid_cells_[layer][x][y];
}

// build_routing_grid Function Implementation
InternalRoutingGrid build_routing_grid(
    const CoreRegion &core_region,
    double lambda_micron,
    double tile_size_lambda,
    int num_layers,
    int default_capacity)
{

    if (lambda_micron <= 0)
    {
        throw std::invalid_argument("Lambda must be positive.");
    }
    if (tile_size_lambda <= 0)
    {
        throw std::invalid_argument("Tile size in lambda must be positive.");
    }
    if (core_region.width_microns <= 0 || core_region.height_microns <= 0)
    {
        throw std::invalid_argument("Core region dimensions must be positive.");
    }
    if (num_layers <= 0)
    {
        throw std::invalid_argument("Number of layers must be positive.");
    }
    if (default_capacity < 0)
    {
        throw std::invalid_argument("Default capacity cannot be negative.");
    }

    double tile_actual_size_microns = tile_size_lambda * lambda_micron;
    if (tile_actual_size_microns == 0)
    { // Avoid division by zero if lambda_micron or tile_size_lambda is extremely small
        throw std::invalid_argument("Calculated tile actual size is zero. Check lambda and tile_size_lambda.");
    }

    int num_x_tiles = static_cast<int>(std::ceil(core_region.width_microns / tile_actual_size_microns));
    int num_y_tiles = static_cast<int>(std::ceil(core_region.height_microns / tile_actual_size_microns));

    if (num_x_tiles <= 0)
        num_x_tiles = 1; // Ensure at least one tile if width is very small but positive
    if (num_y_tiles <= 0)
        num_y_tiles = 1; // Ensure at least one tile if height is very small but positive

    InternalRoutingGrid grid(num_x_tiles, num_y_tiles, num_layers);

    for (int l = 0; l < num_layers; ++l)
    {
        for (int i = 0; i < num_x_tiles; ++i)
        {
            for (int j = 0; j < num_y_tiles; ++j)
            {
                // GridEdge &cell_data = grid.get_grid_cell_data(l, i, j); // Unused variable
                // Or directly: grid.grid_cells_[l][i][j]; if public access is intended here.
                // Let's use direct access since it is public for now for this init function.
                grid.grid_cells_[l][i][j].layer = l;
                grid.grid_cells_[l][i][j].x = i;
                grid.grid_cells_[l][i][j].y = j;
                grid.grid_cells_[l][i][j].demand = 0;                  // Initial demand is zero
                grid.grid_cells_[l][i][j].capacity = default_capacity; // Set to default capacity
            }
        }
    }

    return grid;
}

// --- GlobalRouter Class Implementation ---

#include <queue>     // For std::queue (Soukup)
#include <map>       // For std::map (came_from, cost_so_far)
#include <set>       // For std::set (visited in Soukup)
#include <algorithm> // For std::reverse, std::sort (for priority queue in Hadlock if not using std::priority_queue)
#include <limits>    // For std::numeric_limits
#include <cmath>     // For std::abs (Manhattan distance)

// Helper to create a Pin from grid coordinates
// This might be better placed in Pin struct or as a free function if used elsewhere
Pin make_pin_from_coords(int x, int y, int l, const std::string &net_name = "", const std::string &cell_name = "")
{
    Pin p;
    p.grid_x = x;
    p.grid_y = y;
    p.layer = l;
    p.net_name = net_name;
    p.cell_name = cell_name;
    return p;
}

GlobalRouter::GlobalRouter(InternalRoutingGrid &grid) : routing_grid_(grid) {}

GlobalRoute GlobalRouter::reconstruct_path(const Pin &start_pin, const Pin &end_pin,
                                           const std::map<Pin, Pin> &came_from)
{
    GlobalRoute path;
    Pin current = end_pin;
    std::vector<Pin> pin_path;

    while (!(current.grid_x == start_pin.grid_x && current.grid_y == start_pin.grid_y && current.layer == start_pin.layer))
    {
        pin_path.push_back(current);
        auto it = came_from.find(current);
        if (it == came_from.end())
        {
            // Path is broken or start_pin was the end_pin itself (should be handled before calling)
            return {}; // Return empty path
        }
        current = it->second;
    }
    pin_path.push_back(start_pin); // Add the start pin
    std::reverse(pin_path.begin(), pin_path.end());

    // Convert Pin path to GridEdge path and update demand
    for (size_t i = 0; i < pin_path.size(); ++i)
    {
        path.edges.push_back({pin_path[i].layer, pin_path[i].grid_x, pin_path[i].grid_y, 0, 0}); // demand/cap in GridEdge not directly used from RoutingTileEdge here
        if (i > 0)
        {
            update_edge_demand(pin_path[i - 1], pin_path[i]);
        }
    }
    return path;
}

void GlobalRouter::update_edge_demand(const Pin &from_pin, const Pin &to_pin, int demand_increase)
{
    // Ensure pins are adjacent (differ by 1 in one coord, same layer for now)
    // This simplified version assumes 2D routing on a single layer or layer change handled by via modelling
    // Demand is updated on the G-Cell being entered (to_pin).

    if (from_pin.layer != to_pin.layer)
    {
        // TODO: Handle via demand if multi-layer routing is fully supported
        // This would involve updating demand on the 'via' itself if vias have capacity,
        // or on the cells at (from_pin.x, from_pin.y, to_pin.layer) if changing layers.
        return;
    }

    // Check if to_pin is within grid boundaries before attempting to get its data.
    // get_grid_cell_data() would throw an exception if out of bounds, this is a pre-check.
    if (to_pin.grid_x >= 0 && to_pin.grid_x < routing_grid_.get_num_x_tiles() &&
        to_pin.grid_y >= 0 && to_pin.grid_y < routing_grid_.get_num_y_tiles() &&
        to_pin.layer >= 0 && to_pin.layer < routing_grid_.get_num_layers())
    {
        GridEdge &target_cell_data = routing_grid_.get_grid_cell_data(to_pin.layer, to_pin.grid_x, to_pin.grid_y);
        target_cell_data.demand += demand_increase;

        // Optional: Check for overflow
        // if (target_cell_data.demand > target_cell_data.capacity) {
        //     // Handle overflow (e.g., log, throw exception, or mark as overflowed)
        // }
    }
    // else: to_pin is out of bounds. This should ideally not be reached if get_neighbors is correct
    // and pins are always valid.
}

std::vector<Pin> GlobalRouter::get_neighbors(const Pin &current_pin)
{
    std::vector<Pin> neighbors;
    int x = current_pin.grid_x;
    int y = current_pin.grid_y;
    int l = current_pin.layer;

    // Potential moves (dx, dy, dl)
    int dx[] = {0, 0, 1, -1}; // N, S, E, W
    int dy[] = {1, -1, 0, 0}; // N, S, E, W
    // Not handling layer changes (vias) in this basic get_neighbors for now
    // To add vias: int dl[] = {0,0,0,0,1,-1} and adjust dx/dy to have 0s for via moves

    for (int i = 0; i < 4; ++i)
    {
        int nx = x + dx[i];
        int ny = y + dy[i];
        // int nl = l + dl[i]; // For vias
        int nl = l; // Current: 2D neighbors only

        // Boundary checks for the target neighbor cell (nx, ny, nl)
        if (nx < 0 || nx >= routing_grid_.get_num_x_tiles() ||
            ny < 0 || ny >= routing_grid_.get_num_y_tiles() ||
            nl < 0 || nl >= routing_grid_.get_num_layers())
        {
            continue;
        }

        // Check capacity of the target neighbor cell (nx, ny, nl)
        // get_grid_cell_data will throw if (nl,nx,ny) is out of bounds,
        // but we've already performed boundary checks above.
        const GridEdge &neighbor_cell_data = routing_grid_.get_grid_cell_data(nl, nx, ny);

        if (neighbor_cell_data.demand < neighbor_cell_data.capacity)
        {
            neighbors.push_back(make_pin_from_coords(nx, ny, nl, current_pin.net_name));
        }
    }
    // TODO: Add via capacity checks if nl != l. This would mean checking capacity
    // of the cell (current_pin.x, current_pin.y, new_layer) if vias are modeled as cell transitions.
    return neighbors;
}

GlobalRoute GlobalRouter::route_soukup(const Pin &start_pin, const Pin &end_pin)
{
    std::queue<Pin> q;
    std::map<Pin, Pin> came_from; // Stores parent pointers: came_from[child] = parent
    std::set<Pin> visited_cells;  // To avoid cycles and redundant exploration of cells

    q.push(start_pin);
    visited_cells.insert(start_pin);
    // came_from[start_pin] is not set, reconstruction stops when start_pin is current

    bool found = false;
    while (!q.empty())
    {
        Pin current = q.front();
        q.pop();

        if (current.grid_x == end_pin.grid_x && current.grid_y == end_pin.grid_y && current.layer == end_pin.layer)
        {
            found = true;
            break; // Found the target
        }

        for (const auto &neighbor : get_neighbors(current))
        {
            if (visited_cells.find(neighbor) == visited_cells.end())
            {
                visited_cells.insert(neighbor);
                came_from[neighbor] = current;
                q.push(neighbor);
            }
        }
    }

    if (found)
    {
        return reconstruct_path(start_pin, end_pin, came_from);
    }
    return {}; // No path found
}

// Heuristic for Hadlock (Manhattan distance)
int manhattan_distance(const Pin &p1, const Pin &p2)
{
    return std::abs(p1.grid_x - p2.grid_x) + std::abs(p1.grid_y - p2.grid_y) + std::abs(p1.layer - p2.layer); // Include layer for 3D
}

GlobalRoute GlobalRouter::route_hadlock(const Pin &start_pin, const Pin &end_pin)
{
    // Priority queue stores pairs of (cost, Pin), ordered by cost.
    // Using std::map to simulate priority queue for simplicity here, or std::set for custom comparators.
    // A true std::priority_queue is better: std::priority_queue<std::pair<int, Pin>, std::vector<std::pair<int, Pin>>, std::greater<std::pair<int, Pin>>> pq;
    std::map<int, std::vector<Pin>> priority_queue_map; // cost -> list of pins
    std::map<Pin, Pin> came_from;                       // child -> parent
    std::map<Pin, int> cost_so_far;                     // cost from start (g value - number of detours)

    priority_queue_map[manhattan_distance(start_pin, end_pin)].push_back(start_pin); // Initial cost is purely heuristic (0 detours)
    cost_so_far[start_pin] = 0;

    bool found = false;
    Pin final_pin_at_target; // Store the exact pin that reached the target

    while (!priority_queue_map.empty())
    {
        auto it_lowest_cost_pins = priority_queue_map.begin();
        // int current_priority_val = it_lowest_cost_pins->first; // Unused variable
        Pin current = it_lowest_cost_pins->second.front();
        it_lowest_cost_pins->second.erase(it_lowest_cost_pins->second.begin());
        if (it_lowest_cost_pins->second.empty())
        {
            priority_queue_map.erase(it_lowest_cost_pins);
        }

        if (current.grid_x == end_pin.grid_x && current.grid_y == end_pin.grid_y && current.layer == end_pin.layer)
        {
            found = true;
            final_pin_at_target = current; // Path found using this pin
            break;
        }

        for (const auto &neighbor : get_neighbors(current))
        {
            // Hadlock's cost: prefer straight lines, penalize detours.
            // A simple detour calculation: if neighbor is not on a straight line from current towards end_pin.
            // Cost = cost_so_far[current] + 1 (for segment) + detour_penalty
            int base_cost = 1; // Cost of one segment
            int detour_penalty = 0;

            // Simplified detour: if not moving closer to target by Manhattan distance
            int dist_current_target = manhattan_distance(current, end_pin);
            int dist_neighbor_target = manhattan_distance(neighbor, end_pin);
            if (dist_neighbor_target >= dist_current_target && !(neighbor.grid_x == end_pin.grid_x && neighbor.grid_y == end_pin.grid_y))
            { // Not strictly better and not target
                // This is a very simple detour heuristic. True Hadlock is more about escape points.
                // A more accurate Hadlock cost: g(n) is the number of detours.
                // A detour occurs if the chosen segment does not lie on any shortest path from S to T.
                // Or, more simply, if it moves away from T or parallel when it could move towards.
                Pin came_from_current = current; // Default to current if no parent (start_pin)
                if (came_from.count(current))
                    came_from_current = came_from[current];

                // Check if (came_from_current -> current -> neighbor) is a turn
                bool is_turn = false;
                if (came_from.count(current))
                {
                    Pin prev = came_from.at(current);
                    if (!((neighbor.grid_x - current.grid_x == current.grid_x - prev.grid_x) && (neighbor.grid_y - current.grid_y == current.grid_y - prev.grid_y)))
                    {
                        is_turn = true;
                    }
                }
                if (is_turn)
                    detour_penalty = 1; // Penalize turns strongly
            }

            int new_cost = cost_so_far[current] + base_cost + detour_penalty;

            if (cost_so_far.find(neighbor) == cost_so_far.end() || new_cost < cost_so_far[neighbor])
            {
                cost_so_far[neighbor] = new_cost;
                int priority = new_cost + manhattan_distance(neighbor, end_pin); // f = g + h
                priority_queue_map[priority].push_back(neighbor);
                came_from[neighbor] = current;
            }
        }
    }

    if (found)
    {
        return reconstruct_path(start_pin, final_pin_at_target, came_from);
    }

    return {}; // No path found
}

std::vector<GlobalRoute> GlobalRouter::route_all_nets(const std::vector<Net *> &nets_to_route, const std::string &algorithm)
{
    std::vector<GlobalRoute> all_routes;
    // For multi-pin nets, decomposition is needed (e.g., MST or Steiner tree approximation).
    // For now, let's assume nets are 2-pin or we process the first 2 pins.
    // A more robust approach would be to iterate through pin pairs from an MST of the net.

    for (const auto &net : nets_to_route)
    {
        if (net->nodes.size() < 2)
        {
            continue; // Not routable or already routed (single pin net)
        }

        // Simple decomposition: connect pins sequentially or build MST
        // For this example, connect pins sequentially or build MST
        // For a true global router, one would typically use a Minimum Spanning Tree (MST)
        // to break down multi-pin nets into a series of 2-pin connections.

        // Let's assume Node positions are grid-aligned or can be mapped to grid cells (Pins)
        // This mapping (Node -> Pin) is crucial and depends on how placement output translates to grid.
        // For now, placeholder: use node name as part of pin info.
        // The actual grid_x, grid_y, layer for a Node's pin must be determined from placement data + grid definition.

        // Placeholder: Iterate through pairs of pins for the net
        // This requires Node objects in Net to have accessible, grid-mapped pin locations.
        // For now, this part needs significant refinement based on how Node positions become Pin coords.
        // Let's assume a simple 2-pin net for demonstration if net->nodes.size() == 2.
        if (net->nodes.size() >= 2)
        {
            // TODO: This is a MAJOR simplification. Need to map Node positions to Pin grid coordinates.
            // Assume node->pos is the center and we snap to the nearest grid cell/point.
            // Assume layer 0 for now.
            // The Pin struct should be populated with actual grid coordinates derived from node positions.
            Pin start_pin;
            // start_pin.grid_x = static_cast<int>(net->nodes[0]->pos.x / routing_grid_.tile_width);
            // start_pin.grid_y = static_cast<int>(net->nodes[0]->pos.y / routing_grid_.tile_height);
            // start_pin.layer = 0; // Default layer
            // start_pin.net_name = net->name;
            // start_pin.cell_name = net->nodes[0]->name;

            Pin end_pin;
            // end_pin.grid_x = static_cast<int>(net->nodes[1]->pos.x / routing_grid_.tile_width);
            // end_pin.grid_y = static_cast<int>(net->nodes[1]->pos.y / routing_grid_.tile_height);
            // end_pin.layer = 0;
            // end_pin.net_name = net->name;
            // end_pin.cell_name = net->nodes[1]->name;

            // The following is a conceptual placeholder. The actual pin generation is missing.
            // We need a clear way to get Pins (grid_x, grid_y, layer) for each Node in the Net.
            // For the purpose of having compilable code, we will skip routing if actual pins cannot be derived here.
            // Properly, one would iterate (N-1) times for an N-pin net using an MST.

            if (false)
            { // Guard to prevent execution until Pin creation is solid
                GlobalRoute route;
                if (algorithm == "soukup")
                {
                    route = route_soukup(start_pin, end_pin);
                }
                else if (algorithm == "hadlock")
                {
                    route = route_hadlock(start_pin, end_pin);
                }
                else
                {
                    // Default or throw error
                    route = route_soukup(start_pin, end_pin);
                }

                if (!route.edges.empty())
                {
                    all_routes.push_back(route);
                }
            }
        }
        // For multi-pin nets, loop through (N-1) pairs from MST.
    }
    return all_routes;
}