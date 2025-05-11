#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <list>
#include <algorithm> // For std::min/max used in NetChannelRepresentation merged constructor
#include <limits>    // For std::numeric_limits used in NetChannelRepresentation merged constructor

#define LAMBDA 0.1

class Node;
struct Net;
// Forward declaration for NetChannelRepresentation, needed by Channel
struct NetChannelRepresentation;

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

// Forward declaration for Net, assuming it's defined elsewhere or will be.
// If Net is part of this project and defined in another header,
// that header should be included instead.
struct Net;

struct GridEdge
{
    int layer;    // The layer index of this grid cell/edge segment
    int x;        // The x-coordinate of this grid cell
    int y;        // The y-coordinate of this grid cell
    int demand;   // Current routing demand through this cell/edge
    int capacity; // Maximum routing capacity of this cell/edge
};

struct GlobalRoute
{
    std::vector<GridEdge> edges; // Represents a path as a sequence of grid cells
};

// Definitions for Channel Routing
enum class PinSide
{
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
    UNKNOWN // For unassigned or error states
};

struct ChannelPin
{
    // int pin_id; // Optional: if pins have unique IDs themselves
    int x_coord;  // x-coordinate along the channel length/width
    PinSide side; // Which side of the channel (TOP, BOTTOM, LEFT, RIGHT)
    // int layer;   // Optional: if pins can be on different layers within the channel access point
    // Net* owner_net_ptr; // Optional: direct pointer back to its Net, might be redundant
};

// Represents a net's specific instance and characteristics within a particular channel.
// This object would be created and populated by the ChannelParser.
struct NetChannelRepresentation
{
    Net *original_net;            // Pointer to the global Net object this representation corresponds to
    std::string net_id_str;       // The ID of the net as read from .chn file (e.g., "n0", "n12")
    std::vector<ChannelPin> pins; // List of pins for this net within this channel

    // Fields to be populated and used by the DetailedRouter:
    int left_most_x;      // Leftmost extent (x-coordinate) of the net in the channel
    int right_most_x;     // Rightmost extent (x-coordinate) of the net in the channel
    int track_assignment; // Assigned track number (-1 if not yet assigned)
    bool is_routed;       // Flag indicating if the net (or its current segment) has been routed

    bool is_primary_segment;
    int current_segment_start_x;
    int current_segment_end_x;

    bool is_merged_net;
    std::vector<const NetChannelRepresentation *> constituent_nets;

    NetChannelRepresentation *parent_representation; // For dogleg segments, points to the original NetChannelRepresentation

    NetChannelRepresentation(Net *net_ptr = nullptr, const std::string &id_str = "") : original_net(net_ptr),
                                                                                       net_id_str(id_str),
                                                                                       left_most_x(-1),
                                                                                       right_most_x(-1),
                                                                                       track_assignment(-1),
                                                                                       is_routed(false),
                                                                                       is_primary_segment(true),
                                                                                       current_segment_start_x(-1),
                                                                                       current_segment_end_x(-1),
                                                                                       is_merged_net(false),
                                                                                       parent_representation(nullptr) // Initialize new member
    {
        if (net_ptr && id_str.empty())
        {
            net_id_str = net_ptr->name;
        }
    }

    NetChannelRepresentation(const std::vector<const NetChannelRepresentation *> &members, const std::string &merged_id) : original_net(nullptr),
                                                                                                                           net_id_str(merged_id),
                                                                                                                           left_most_x(std::numeric_limits<int>::max()),
                                                                                                                           right_most_x(std::numeric_limits<int>::min()),
                                                                                                                           track_assignment(-1),
                                                                                                                           is_routed(false),
                                                                                                                           is_primary_segment(true),
                                                                                                                           is_merged_net(true),
                                                                                                                           constituent_nets(members),
                                                                                                                           parent_representation(nullptr) // Initialize new member
    {
        if (members.empty())
        {
            current_segment_start_x = left_most_x;
            current_segment_end_x = right_most_x;
            return;
        }
        for (const auto *member_net : members)
        {
            if (member_net)
            {
                pins.insert(pins.end(), member_net->pins.begin(), member_net->pins.end());
                left_most_x = std::min(left_most_x, member_net->left_most_x);
                right_most_x = std::max(right_most_x, member_net->right_most_x);
            }
        }
        if (left_most_x == std::numeric_limits<int>::max())
        {
            left_most_x = -1;
            right_most_x = -1;
        }
        current_segment_start_x = left_most_x;
        current_segment_end_x = right_most_x;
        if (!members.empty() && members[0] != nullptr)
        {
            original_net = members[0]->original_net;
        }
    }

    // New constructor for dogleg segments
    NetChannelRepresentation(NetChannelRepresentation *parent, const std::string &segment_id_str, int start_x, int end_x, const std::vector<ChannelPin> &segment_pins) : original_net(parent ? parent->original_net : nullptr),
                                                                                                                                                                         net_id_str(segment_id_str),
                                                                                                                                                                         pins(segment_pins),
                                                                                                                                                                         left_most_x(start_x),
                                                                                                                                                                         right_most_x(end_x),
                                                                                                                                                                         track_assignment(-1),
                                                                                                                                                                         is_routed(false),
                                                                                                                                                                         is_primary_segment(false), // Dogleg segments are not primary
                                                                                                                                                                         current_segment_start_x(start_x),
                                                                                                                                                                         current_segment_end_x(end_x),
                                                                                                                                                                         is_merged_net(parent ? parent->is_merged_net : false), // Inherit merged status
                                                                                                                                                                         // constituent_nets: if parent was merged, this segment is part of that merged concept. Could copy or leave empty.
                                                                                                                                                                         // For simplicity, dogleg segments don't re-list constituents of a merged parent directly.
                                                                                                                                                                         parent_representation(parent)
    {
        // If the parent was a merged net, the dogleg segment is still part of that conceptual merged net.
        if (parent && parent->is_merged_net)
        {
            constituent_nets = parent->constituent_nets; // Share constituent info
        }
    }
};

struct Channel
{
    int id;                                       // Unique identifier for the channel
    int width;                                    // Width of the channel (e.g., in lambda units or grid units - channel length)
    int height;                                   // Height of the channel (e.g., number of tracks or physical height)
    std::vector<NetChannelRepresentation *> nets; // Pointers to channel-specific net representations
    // Note: Storing raw pointers implies manual memory management for NetChannelRepresentation objects
    // or that these are non-owning pointers to objects managed elsewhere (e.g., by ChannelParser or a global store).
};

#endif // STRUCTURES_H
