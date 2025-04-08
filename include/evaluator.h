#pragma once

#include "structures.h"
#include <vector>
#include <unordered_map>
#include <string>

// Forward declaration
class Circuit;

class Evaluator
{
public:
    // Constructor takes the final circuit and placement results
    Evaluator(Circuit &circuit, const std::unordered_map<std::string, Point> &placement);

    // --- Core Evaluation Metrics ---

    // Calculate total Half-Perimeter Wirelength (HPWL)
    double calculate_total_hpwl();

    // Calculate the bounding box area of the placed movable nodes
    double calculate_placement_area();

    // Estimate total routing overflow based on the grid and estimated demand
    // Returns a single metric (e.g., sum of overflows over capacity across all edges)
    double estimate_total_overflow();

    // Estimate critical path delay (requires timing information in Circuit)
    // Placeholder - Requires a proper timing model (e.g., Elmore delay)
    double estimate_critical_delay();

    // --- Feedback Generation ---

    // Identify regions with high estimated routing congestion
    // Returns a list of descriptive strings or region identifiers
    std::vector<std::string> get_congested_regions_feedback();

private:
    Circuit &circuit_ref;                                        // Reference to the circuit data
    const std::unordered_map<std::string, Point> &placement_ref; // Reference to final placement

    // Helper to calculate HPWL for a single net
    double calculate_net_hpwl(const Net &net);

    // Helper to estimate routing demand on the grid based on net bounding boxes
    // This is a crucial but complex step, often involving approximations.
    void estimate_routing_demand();

    // Helper for delay calculation (e.g., implementing Elmore delay)
    // Needs access to pin offsets, net topology, and potentially wire resistance/capacitance
    double calculate_elmore_delay(const Net &net);

    // Internal state for calculated metrics
    bool demand_estimated = false;                        // Flag to avoid re-calculating demand
    std::map<std::pair<int, int>, double> routing_demand; // Grid cell -> demand mapping
};