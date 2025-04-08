#pragma once

#include "structures.h"
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <limits>

// Forward declaration
class Circuit;

class FMPartitioner
{
public:
    // Constructor takes a reference to the main circuit data
    FMPartitioner(Circuit &circuit);

    // Main partitioning function
    // Takes a list of nodes to partition (pointers)
    // Returns the two partitions and the cut size
    PartitionResult partition(const std::vector<Node *> &nodes_to_partition,
                              double balance_factor_min, // Min area ratio for partition A (e.g., 0.45)
                              double balance_factor_max);

private:
    Circuit &circuit_ref; // Reference to the overall circuit data

    // Internal FM data structures
    // Gain bucket: Vector index = gain + p_max
    std::vector<std::list<Node *>> gain_bucket;
    int min_gain_index = std::numeric_limits<int>::max();
    int max_gain_index = -std::numeric_limits<int>::max();
    int p_max = 0; // Max degree of a node, needed for gain calculation range

    // State for the current pass
    double current_target_area_A = 0.0;
    double current_target_area_B = 0.0;
    double current_area_A = 0.0;
    double current_area_B = 0.0;
    double balance_factor_min_d = 0.0;
    double balance_factor_max_d = 0.0;

    std::vector<Node *> current_nodes; // Nodes being partitioned in this call

    // Constants
    static const int MAX_ITERATIONS = 10;
    static const int MIN_NODES_PER_REGION = 5;

    // --- Internal FM Methods ---

    // 1. Initialization
    void initialize_partition(const std::vector<Node *> &nodes, double target_balance);

    void initialize_gains_and_buckets();
    void printGainBucket();
    void reset_fm_state();

    // 2. Main Pass Logic
    PartitionResult run_fm_pass();

    bool check_area_balance(Node *node_to_move);
    // void move_node(Node *node_to_move);

    // 3. Gain Update
    void update_gains_after_move(Node *moved_node);
    void update_single_node_gain(Node *node, int delta_gain); // Helper to update gain & bucket

    // 4. Utility
    int calculate_cut_size();

    double calculate_total_area(const std::vector<Node *> &nodes);

    // New methods for placement objectives
    double calculate_total_wire_length();
    double estimate_congestion();
    double estimate_timing();
    bool is_better_partition(int current_cut, double current_wire_length,
                             double current_congestion, double current_timing,
                             int best_cut, double best_wire_length,
                             double best_congestion, double best_timing);

    // Rollback method
    void rollback_to_best_cut(const std::vector<Node *> &move_sequence,
                              const std::vector<int> &cut_sizes,
                              const std::vector<double> &area_A_history,
                              const std::vector<double> &area_B_history,
                              int best_move_index, PartitionResult best_pass_result);
};