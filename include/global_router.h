#ifndef GLOBAL_ROUTER_H
#define GLOBAL_ROUTER_H

#include <vector>
#include <string>
#include "structures.h" // Contains RoutingGrid, GlobalRoute, GridEdge, Net

// Forward declarations (if not fully defined in structures.h or if prefered)
// struct RoutingGrid; // Already in structures.h
// struct Net;         // Already in structures.h
// struct GlobalRoute; // Already in structures.h
// struct GridEdge;    // Already in structures.h
struct Pin; // Assuming a Pin structure to define start/end points of a route

// Represents a pin on a cell, used as a target for routing.
// A pin has a location on the grid.
struct Pin
{
    std::string net_name;
    std::string cell_name; // Optional: name of the cell this pin belongs to
    int grid_x;            // x-coordinate in the routing grid
    int grid_y;            // y-coordinate in the routing grid
    int layer;             // layer of the pin (important for 3D routing)

    // Helper to compare pins, e.g. for map keys or set elements
    bool operator<(const Pin &other) const
    {
        if (grid_x != other.grid_x)
            return grid_x < other.grid_x;
        if (grid_y != other.grid_y)
            return grid_y < other.grid_y;
        return layer < other.layer;
    }
};

// Placeholder for CoreRegion structure.
// This would typically define the bounding box and origin of the routing area.
struct CoreRegion
{
    double x_min_microns;  // Minimum x-coordinate of the core region
    double y_min_microns;  // Minimum y-coordinate of the core region
    double width_microns;  // Width of the core routing area in microns
    double height_microns; // Height of the core routing area in microns
};

class InternalRoutingGrid
{
public:
    // Constructor
    InternalRoutingGrid(int num_x_tiles, int num_y_tiles, int num_layers);

    // Accessors
    int get_num_x_tiles() const { return num_x_tiles_; }
    int get_num_y_tiles() const { return num_y_tiles_; }
    int get_num_layers() const { return num_layers_; }

    // Get a modifiable reference to a specific GridEdge (G-Cell)
    GridEdge &get_grid_cell_data(int layer, int x, int y);
    // Get a const reference to a specific GridEdge (G-Cell)
    const GridEdge &get_grid_cell_data(int layer, int x, int y) const;

    // Potentially methods to get neighbors, capacities, etc. can be added later.

    // Making grid_cells_ public for easier direct manipulation in build_routing_grid initially,
    // but ideally, it would be private with more specific public methods.
    // For the purpose of initialization and directness per plan:
    // "std::vector<GridEdge> edges;" in GlobalRoute and GridEdge struct itself,
    // suggests GridEdge is the primary data unit.
    // So, RoutingGrid will own a 3D vector of these.
    std::vector<std::vector<std::vector<GridEdge>>> grid_cells_;

private:
    int num_x_tiles_;
    int num_y_tiles_;
    int num_layers_;
};

// Function to build the routing grid
InternalRoutingGrid build_routing_grid(
    const CoreRegion &core_region, // Defines the total area for routing
    double lambda_micron,          // The size of one lambda unit in microns
    double tile_size_lambda,       // The size of one grid tile in lambda units (e.g., 10 means tile is 10λ x 10λ)
    int num_layers,                // Number of routing layers
    int default_capacity           // Default capacity for each G-Cell's edges
);

class GlobalRouter
{
public:
    GlobalRouter(InternalRoutingGrid &grid); // Constructor taking the routing grid

    // Performs Soukup's algorithm (BFS-based maze routing) for a single 2-pin connection.
    // Updates the demand on grid edges.
    GlobalRoute route_soukup(const Pin &start_pin, const Pin &end_pin);

    // Performs Hadlock's algorithm (A*-based, cost = detours) for a single 2-pin connection.
    // Updates the demand on grid edges.
    GlobalRoute route_hadlock(const Pin &start_pin, const Pin &end_pin);

    // Routes all nets in a given list.
    // Each net will be decomposed into 2-pin connections (e.g., using a minimum spanning tree).
    // This is a higher-level function that will call a 2-pin routing algorithm for each segment.
    std::vector<GlobalRoute> route_all_nets(const std::vector<Net *> &nets_to_route, const std::string &algorithm = "soukup");

private:
    InternalRoutingGrid &routing_grid_; // Reference to the main routing grid

    // Helper function to reconstruct the path from parent pointers (used by BFS/A*)
    GlobalRoute reconstruct_path(const Pin &start_pin, const Pin &end_pin,
                                 const std::map<Pin, Pin> &came_from);

    // Helper to get neighboring grid cells/edges for pathfinding
    std::vector<Pin> get_neighbors(const Pin &current_pin);

    // Helper to update demand on a specific edge in the routing_grid_
    // The GridEdge here would represent the 'to' node of an edge from a 'from' node.
    void update_edge_demand(const Pin &from_pin, const Pin &to_pin, int demand_increase = 1);
};

#endif // GLOBAL_ROUTER_H