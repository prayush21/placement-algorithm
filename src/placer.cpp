#include "../include/placer.h"
#include <iostream>
#include <algorithm>

// Constructor implementation
Placer::Placer(Circuit &circuit, FMPartitioner &fm)
    : circuit_ref(circuit), fm_ref(fm)
{
    // Initialize any necessary state
}

// Main placement function
void Placer::place(PlacementStrategy strategy)
{
    std::cout << "[Placer] Starting placement with strategy: " << strategy << std::endl;

    // Create initial region covering the entire core area
    PlacementRegion initial_region(circuit_ref.core_region);

    // Collect movable nodes
    std::vector<Node *> movable_nodes;
    for (auto &node_pair : circuit_ref.cell_map)
    {
        if (node_pair.second.type == MOVABLE)
        {
            movable_nodes.push_back(&node_pair.second);
        }
    }
    initial_region.assigned_nodes = movable_nodes;

    // Call appropriate recursive placement function based on strategy
    switch (strategy)
    {
    case BISECTION:
        recursive_bisection(initial_region, 0);
        break;
    case QUADRATURE:
        recursive_quadrature(initial_region, 0);
        break;
    case SLICE_BISECTION:
        recursive_slice_bisection(initial_region, 0, true);
        break;
    case CUT_ORIENTED:
        recursive_cut_oriented(initial_region, 0);
        break;
    default:
        throw std::runtime_error("Unknown placement strategy");
    }
}

const std::unordered_map<std::string, Point> &Placer::get_placement() const
{
    return final_placement;
}

// --- Private Method Implementations ---

void Placer::recursive_bisection(PlacementRegion current_region, int level)
{
    bool should_log = level + 1;
    if (should_log)
    {
        std::cout << "[Placer][Level " << level << "] Quadrature on region with " << current_region.assigned_nodes.size() << " nodes." << std::endl;
        std::cout << "Approximately " << nodes_placed << " Nodes placed.\n";
    }
    if (should_terminate(current_region, level))
    {
        assign_coords_in_leaf(current_region);
        return;
    }

    // Define balance constraints for FM partitioning
    const double balance_min = 0.45;
    const double balance_max = 0.55;

    // Determine cut direction (alternate: horizontal for even levels, vertical for odd)
    bool cut_vertically = (level % 2 != 0);

    // Partition nodes
    PartitionResult split = fm_ref.partition(current_region.assigned_nodes, balance_min, balance_max);

    // Create sub-regions
    Rect region1_bounds, region2_bounds;
    double center_x = (current_region.bounds.bottom_left.x + current_region.bounds.top_right.x) / 2.0;
    double center_y = (current_region.bounds.bottom_left.y + current_region.bounds.top_right.y) / 2.0;

    if (cut_vertically)
    {
        // Vertical cut
        region1_bounds = {current_region.bounds.bottom_left, {center_x, current_region.bounds.top_right.y}}; // Left
        region2_bounds = {{center_x, current_region.bounds.bottom_left.y}, current_region.bounds.top_right}; // Right
    }
    else
    {
        // Horizontal cut
        region1_bounds = {current_region.bounds.bottom_left, {current_region.bounds.top_right.x, center_y}}; // Bottom
        region2_bounds = {{current_region.bounds.bottom_left.x, center_y}, current_region.bounds.top_right}; // Top
    }

    PlacementRegion region1(region1_bounds);
    region1.assigned_nodes = split.partition_A;

    PlacementRegion region2(region2_bounds);
    region2.assigned_nodes = split.partition_B;

    // Recursive calls
    recursive_bisection(region1, level + 1);
    recursive_bisection(region2, level + 1);
}

void Placer::recursive_quadrature(PlacementRegion current_region, int level)
{
    // Print progress every 6 levels
    bool should_log = level + 1;

    if (should_log)
    {
        std::cout << "[Placer][Level " << level << "] Quadrature on region with " << current_region.assigned_nodes.size() << " nodes." << std::endl;
        std::cout << "Approximately " << nodes_placed << " Nodes placed.\n";
    }

    if (should_terminate(current_region, level))
    {
        assign_coords_in_leaf(current_region);
        return;
    }

    // Define balance constraints for FM partitioning
    const double balance_min = 0.45;
    const double balance_max = 0.55;

    // Calculate center coordinates for geometric splitting
    double center_x = (current_region.bounds.bottom_left.x + current_region.bounds.top_right.x) / 2.0;
    double center_y = (current_region.bounds.bottom_left.y + current_region.bounds.top_right.y) / 2.0;

    // --- Hierarchical Partitioning Steps ---

    PartitionResult first_split = fm_ref.partition(current_region.assigned_nodes, balance_min, balance_max);
    // std::cout << "Print Region 1: ";
    // for (auto node : first_split.partition_A)
    // {
    //     std::cout << node->name << " ,";
    // }
    // std::cout << std::endl;
    // std::cout << "Print Region 2: ";
    // for (auto node : first_split.partition_B)
    // {
    //     std::cout << node->name << " ,";
    // }
    // std::cout << std::endl;
    std::vector<Node *> group1_nodes = first_split.partition_A; // Let's call this group 1
    std::vector<Node *> group2_nodes = first_split.partition_B; // Let's call this group 2

    PartitionResult group1_split = fm_ref.partition(group1_nodes, balance_min, balance_max);
    std::vector<Node *> bottom_left_nodes = group1_split.partition_A;
    std::vector<Node *> top_left_nodes = group1_split.partition_B;

    PartitionResult group2_split = fm_ref.partition(group2_nodes, balance_min, balance_max);
    std::vector<Node *> bottom_right_nodes = group2_split.partition_A;
    std::vector<Node *> top_right_nodes = group2_split.partition_B;

    // --- Define Geometric Sub-Regions (Quadrants) ---
    Rect bl_bounds = {current_region.bounds.bottom_left, {center_x, center_y}};
    Rect tl_bounds = {{current_region.bounds.bottom_left.x, center_y}, {center_x, current_region.bounds.top_right.y}};
    Rect br_bounds = {{center_x, current_region.bounds.bottom_left.y}, {current_region.bounds.top_right.x, center_y}};
    Rect tr_bounds = {{center_x, center_y}, current_region.bounds.top_right};

    // --- Create Sub-Regions and Assign Node Groups ---
    PlacementRegion bottom_left_region(bl_bounds);
    bottom_left_region.assigned_nodes = bottom_left_nodes;

    PlacementRegion top_left_region(tl_bounds);
    top_left_region.assigned_nodes = top_left_nodes;

    PlacementRegion bottom_right_region(br_bounds);
    bottom_right_region.assigned_nodes = bottom_right_nodes;

    PlacementRegion top_right_region(tr_bounds);
    top_right_region.assigned_nodes = top_right_nodes;

    // --- Recursive Calls ---
    recursive_quadrature(bottom_left_region, level + 1);
    recursive_quadrature(top_left_region, level + 1);
    recursive_quadrature(bottom_right_region, level + 1);
    recursive_quadrature(top_right_region, level + 1);
}

void Placer::recursive_slice_bisection(PlacementRegion current_region, int level, bool vertical_first)
{
    bool should_log = level + 1;
    if (should_log)
    {
        std::cout << "[Placer][Level " << level << "] Quadrature on region with " << current_region.assigned_nodes.size() << " nodes." << std::endl;
        std::cout << "Approximately " << nodes_placed << " Nodes placed.\n";
    }

    bool cut_vertically = vertical_first;

    if (should_terminate(current_region, level))
    {
        assign_coords_in_leaf(current_region);
        return;
    }

    // Define balance constraints
    const double balance_min = 0.45;
    const double balance_max = 0.55;

    // Partition nodes
    PartitionResult split = fm_ref.partition(current_region.assigned_nodes, balance_min, balance_max);

    // Create sub-regions based on cut direction
    Rect region1_bounds, region2_bounds;
    double center_x = (current_region.bounds.bottom_left.x + current_region.bounds.top_right.x) / 2.0;
    double center_y = (current_region.bounds.bottom_left.y + current_region.bounds.top_right.y) / 2.0;

    if (cut_vertically)
    {
        // Vertical cut
        region1_bounds = {current_region.bounds.bottom_left, {center_x, current_region.bounds.top_right.y}}; // Left
        region2_bounds = {{center_x, current_region.bounds.bottom_left.y}, current_region.bounds.top_right}; // Right
    }
    else
    {
        // Horizontal cut
        region1_bounds = {current_region.bounds.bottom_left, {current_region.bounds.top_right.x, center_y}}; // Bottom
        region2_bounds = {{current_region.bounds.bottom_left.x, center_y}, current_region.bounds.top_right}; // Top
    }

    PlacementRegion region1(region1_bounds);
    region1.assigned_nodes = split.partition_A;

    PlacementRegion region2(region2_bounds);
    region2.assigned_nodes = split.partition_B;

    // Recursive calls: Alternate the cut direction for the next level
    recursive_slice_bisection(region1, level + 1, !cut_vertically);
    recursive_slice_bisection(region2, level + 1, !cut_vertically);
}

void Placer::recursive_cut_oriented(PlacementRegion current_region, int level)
{
    bool should_log = level + 1;
    if (should_log)
    {
        std::cout << "[Placer][Level " << level << "] Quadrature on region with " << current_region.assigned_nodes.size() << " nodes." << std::endl;
        std::cout << "Approximately " << nodes_placed << " Nodes placed.\n";
    }
    if (should_terminate(current_region, level))
    {
        assign_coords_in_leaf(current_region);
        return;
    }

    // Define balance constraints
    const double balance_min = 0.45;
    const double balance_max = 0.55;

    // Determine cut direction based on aspect ratio
    double width = current_region.bounds.top_right.x - current_region.bounds.bottom_left.x;
    double height = current_region.bounds.top_right.y - current_region.bounds.bottom_left.y;
    bool cut_vertically = (width >= height); // Cut vertically if width is greater or equal

    // Partition nodes
    PartitionResult split = fm_ref.partition(current_region.assigned_nodes, balance_min, balance_max);
    // Create sub-regions
    Rect region1_bounds, region2_bounds;
    double center_x = (current_region.bounds.bottom_left.x + current_region.bounds.top_right.x) / 2.0;
    double center_y = (current_region.bounds.bottom_left.y + current_region.bounds.top_right.y) / 2.0;

    if (cut_vertically)
    {
        // Vertical cut
        region1_bounds = {current_region.bounds.bottom_left, {center_x, current_region.bounds.top_right.y}}; // Left
        region2_bounds = {{center_x, current_region.bounds.bottom_left.y}, current_region.bounds.top_right}; // Right
    }
    else
    {
        // Horizontal cut
        region1_bounds = {current_region.bounds.bottom_left, {current_region.bounds.top_right.x, center_y}}; // Bottom
        region2_bounds = {{current_region.bounds.bottom_left.x, center_y}, current_region.bounds.top_right}; // Top
    }

    PlacementRegion region1(region1_bounds);
    region1.assigned_nodes = split.partition_A;

    PlacementRegion region2(region2_bounds);
    region2.assigned_nodes = split.partition_B;

    // Recursive calls
    recursive_cut_oriented(region1, level + 1);
    recursive_cut_oriented(region2, level + 1);
}

void Placer::assign_coords_in_leaf(PlacementRegion &leaf_region)
{
    // Calculate the center of the leaf region
    double center_x = (leaf_region.bounds.bottom_left.x + leaf_region.bounds.top_right.x) / 2.0;
    double center_y = (leaf_region.bounds.bottom_left.y + leaf_region.bounds.top_right.y) / 2.0;
    Point center_point = {center_x, center_y};

    // Assign the center point to all nodes in this leaf region
    for (Node *node : leaf_region.assigned_nodes)
    {
        if (node != nullptr)
        {
            // Assign the calculated center point
            nodes_placed++;
            final_placement[node->name] = center_point;
        }
    }
}

bool Placer::should_terminate(const PlacementRegion &region, int level)
{
    // Check recursion depth
    if (level >= MAX_RECURSION_DEPTH)
    {
        return true;
    }

    // Check number of nodes
    if (region.assigned_nodes.size() <= MIN_NODES_PER_REGION)
    {
        return true;
    }

    // Check region size relative to core area
    double region_area = (region.bounds.top_right.x - region.bounds.bottom_left.x) *
                         (region.bounds.top_right.y - region.bounds.bottom_left.y);
    double core_area = (circuit_ref.core_region.top_right.x - circuit_ref.core_region.bottom_left.x) *
                       (circuit_ref.core_region.top_right.y - circuit_ref.core_region.bottom_left.y);

    // if (region_area / core_area < MIN_REGION_AREA_RATIO)
    // {
    //     return true;
    // }

    return false;
}