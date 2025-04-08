#ifndef PLACER_H
#define PLACER_H

#include "structures.h"
#include "fm.h" // Needs the partitioner
#include <vector>
#include <map>
#include <string>

// Forward declarations
class Circuit;
class FMPartitioner;

// Enum for different placement strategies
enum PlacementStrategy
{
    BISECTION,       // Basic recursive bisection
    QUADRATURE,      // Divide into four quadrants
    SLICE_BISECTION, // Vertical then horizontal slices (or vice versa)
    CUT_ORIENTED     // Cut orientation depends on region aspect ratio
    // Add more strategies as needed
};

// Represents a region during recursive placement
struct PlacementRegion
{
    Rect bounds;                        // Geometric boundaries of the region
    double target_capacity;             // Target area capacity (sum of node areas assigned)
    double current_area = 0.0;          // Actual area occupied by nodes assigned so far
    std::vector<Node *> assigned_nodes; // Nodes belonging to this region

    // Constructor
    PlacementRegion(Rect b = {{0, 0}, {0, 0}}, double cap = 0.0) : bounds(b), target_capacity(cap) {}
};

// Type alias for terminal anchor information
// Maps an internal node to a list of external anchor points it connects to.
using TerminalAnchors = std::map<Node *, std::vector<Point>>;

class Placer
{
public:
    // Constructor requires the circuit and a partitioner instance
    Placer(Circuit &circuit, FMPartitioner &fm);

    // Main placement function
    void place(PlacementStrategy strategy);

    // Get the final calculated placement coordinates
    const std::map<std::string, Point> &get_placement() const;

private:
    Circuit &circuit_ref;                         // Reference to the circuit data
    FMPartitioner &fm_ref;                        // Reference to the FM partitioner
    std::map<std::string, Point> final_placement; // Node name -> bottom-left coord

    // Recursive placement functions for each strategy
    void recursive_bisection(PlacementRegion current_region, int level);
    void recursive_quadrature(PlacementRegion current_region, int level);
    void recursive_slice_bisection(PlacementRegion current_region, int level, bool vertical_first);
    void recursive_cut_oriented(PlacementRegion current_region, int level);
    // ... add functions for other strategies

    // Helper to assign final coordinates within leaf regions (simple packing/spreading)
    void assign_coords_in_leaf(PlacementRegion &leaf_region);

    // Termination condition checker for recursion
    bool should_terminate(const PlacementRegion &region, int level);

    // Helper for terminal propagation
    TerminalAnchors calculate_terminal_anchors(const PlacementRegion &current_region);

    // Constants for termination (examples)
    const int MAX_RECURSION_DEPTH = 10;
    const int MIN_NODES_PER_REGION = 1;
    // const double MIN_REGION_AREA_RATIO = 0.01; // Compared to total core area
};

#endif // PLACER_H
