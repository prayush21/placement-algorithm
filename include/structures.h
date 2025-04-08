#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <list>

#define LAMBDA 0.1

class Node;
class Net;

struct Point
{
    double x = 0.0;
    double y = 0.0;
};

struct Rect
{
    Point bottom_left;
    Point top_right;
};

enum NodeType
{
    MOVABLE,
    FIXED,
    TERMINAL,
    TERMINAL_NI
};

struct Node
{
    std::string name;
    NodeType type = MOVABLE; // Default type
    double width = 0.0;
    double height = 0.0;
    double area = 0.0;
    Point pos;             // Final or fixed position
    int partition_id = -1; // Temp assignment during FM, -1 means unassigned

    // For non-rectangular fixed nodes:
    std::vector<Rect> shapes; // If complex shape

    // FM-related data (inspired by fm_partitioning.cpp)
    int fm_gain = 0;
    bool fm_locked = false;
    // Iterator for gain bucket optimization (if using list-based bucket)
    std::list<Node *>::iterator bucket_iterator; // Added: Iterator pointing to this node in the gain bucket list

    // List of connected nets (could store indices or pointers)
    std::vector<Net *> nets; // Pointers to connected nets

    // Constructor
    Node(const std::string &n = "", NodeType t = MOVABLE, double w = 0.0, double h = 0.0)
        : name(n), type(t), width(w), height(h), area(w * h) {}

    Node(const std::string &n, double w, double h)
        : name(n), width(w), height(h) {}

    Node(const std::string &n, double w, double h, double area)
        : name(n), width(w), height(h), area(area) {}
};

struct Net
{
    std::string name;
    std::vector<Node *> nodes; // Pointers to nodes in this net
    double weight = 1.0;

    bool is_critical = false;

    // FM-related data (inspired by fm_partitioning.cpp)
    int partition_A_count = 0; // Count of nodes in partition A
    int partition_B_count = 0; // Count of nodes in partition B

    // Constructor
    Net(const std::string &n = "") : name(n) {}
};

struct RoutingTileEdge
{
    int capacity = 0;
    int demand = 0;                      // Estimated after placement
    double wire_width_spacing_sum = 0.0; // From .route MinWireWidth+MinWireSpacing
};

struct RoutingGrid
{
    int grid_x = 0;
    int grid_y = 0;
    int num_layers = 0;
    Point origin;
    double tile_width = LAMBDA;
    double tile_height = LAMBDA;

    // Store edge capacities/demand (complex: 3D vector or map)
    // Map layer -> Map edge_coord -> RoutingTileEdge
    // Using vectors for simplicity, assuming dense grid
    std::vector<std::vector<std::vector<RoutingTileEdge>>> horizontal_edges; // [layer][y][x]
    std::vector<std::vector<std::vector<RoutingTileEdge>>> vertical_edges;   // [layer][y][x]

    // Store blockages (could be represented in edge capacities or a separate structure)
    // std::vector<Rect> blockages; // Example
};

struct CoreRow
{
    int id;
    double y_coord; // Bottom y-coordinate of the row
    double x_coord; // Left x-coordinate of the row
    double num_sites;
    double site_width;
    double site_spacing;

    double height;
    double current_width; // Accumulated width of cells placed in this row
    std::vector<std::string> cells_in_row;

    // Constructor
    CoreRow(int id = 0, double y_coord = 0.0, double x_coord = 0.0, double num_sites = 0.0, double site_width = 0.0, double site_spacing = 0.0, double height = 0.0, double current_width = 0.0)
        : id(id), y_coord(y_coord), x_coord(x_coord), num_sites(num_sites), site_width(site_width), site_spacing(site_spacing), height(height), current_width(current_width) {}
};

// struct Region
// {
//     double x_min, y_min, x_max, y_max;
//     double width() const { return x_max - x_min; }
//     double height() const { return y_max - y_min; }
//     double area() const { return width() * height(); }
//     // Constructor
//     Region(double x0 = 0, double y0 = 0, double x1 = 0, double y1 = 0) : x_min(x0), y_min(y0), x_max(x1), y_max(y1) {}
// };

struct Circuit
{
    std::unordered_map<std::string, Node> cell_map; // Renamed from nodes to cell_map
    std::unordered_map<std::string, Net> net_map;   // Renamed from nets to net_map
    Rect core_region;
    std::vector<CoreRow> core_rows;
    RoutingGrid routing_grid;
    // Other info from .aux, .scl etc.
    std::string benchmark_name;

    // Store technology info (.tech LEF/DEF info if parsed)

    // Helper to add a node safely
    Node *add_node(const std::string &name, NodeType type, double w, double h);
    // Helper to add a net safely
    Net *add_net(const std::string &name);
};

// Result structure for the partitioner
struct PartitionResult
{
    std::vector<Node *> partition_A; // Pointers to nodes in partition A
    std::vector<Node *> partition_B; // Pointers to nodes in partition B
    int cut_size = -1;               // -1 indicates invalid/not computed
};

#endif // STRUCTURES_H
