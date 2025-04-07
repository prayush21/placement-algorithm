#include "../include/fm.h"
#include "../include/structures.h"
#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <limits>
#include <cmath>
#include <cstdlib> // for rand(), srand()
#include <ctime>   // for time()
#include <fstream>

// Constructor implementation
FMPartitioner::FMPartitioner(Circuit &circuit) : circuit_ref(circuit)
{
    // Initialize any necessary state
}

// Main partitioning function
PartitionResult FMPartitioner::partition(const std::vector<Node *> &nodes_to_partition,
                                         double balance_factor_min,
                                         double balance_factor_max)
{
    // std::cout << "[FM] Starting partitioning with " << nodes_to_partition.size() << " nodes" << std::endl;

    // Store nodes to partition
    current_nodes = nodes_to_partition;
    if (current_nodes.empty())
    {
        // std::cout << "[FM] No nodes provided for partitioning. Returning empty result." << std::endl;
        return PartitionResult{}; // Return default/empty result
    }

    // Calculate total area and set target areas based on the nodes being partitioned
    double total_area_partitioned = calculate_total_area(current_nodes);
    // std::cout << "Total Area to Partition" << total_area_partitioned << std::endl;
    // Assuming the balance factor applies to the area of the nodes being partitioned
    current_target_area_A = total_area_partitioned * balance_factor_min;
    current_target_area_B = total_area_partitioned * balance_factor_max;

    // std::cout << "[FM] Target Area Range: [" << current_target_area_A << ", " << current_target_area_B << "]" << std::endl;

    // Initialize the partition randomly for the nodes to be partitioned
    initialize_partition(current_nodes, (balance_factor_min + balance_factor_max) / 2.0); // Use average balance for initial random split

    // Store the best overall partition found
    PartitionResult best_result_overall;
    best_result_overall.cut_size = std::numeric_limits<int>::max();
    std::unordered_map<Node *, int> best_partition_map;

    int max_passes = 10; // Maximum number of passes
    bool improvement_found = true;

    for (int pass = 1; pass <= max_passes && improvement_found; ++pass)
    {
        // std::cout << "\n[FM] --- Starting Pass " << pass << " ---" << std::endl;
        improvement_found = false; // Assume no improvement in this pass initially

        // Calculate initial cut size for this pass
        int cut_before_pass = calculate_cut_size();

        // std::cout << "[FM] Cut size before pass " << pass << ": " << cut_before_pass << std::endl;

        // Store current partition state in case this pass is the best
        if (cut_before_pass < best_result_overall.cut_size)
        {
            // std::cout << "[FM] Storing current partition as potentially best (Cut: " << cut_before_pass << ")." << std::endl;
            best_result_overall.cut_size = cut_before_pass;
            best_partition_map.clear();
            for (Node *node : current_nodes)
            {
                best_partition_map[node] = node->partition_id;
            }
            // Also store current areas A and B if needed in result
            // best_result_overall.area_A = current_area_A;
            // best_result_overall.area_B = current_area_B;
            // std::cout << "Current(Best) Area A:" << current_area_A << " :";
            // for (Node *node : current_nodes)
            // {
            //     if (node->partition_id == 0)
            //         std::cout << node->name << ", ";
            // }
            // std::cout << std::endl;
            // std::cout << "Current(Best) Area B:" << current_area_B << " :";

            // for (Node *node : current_nodes)
            // {
            //     if (node->partition_id == 1)
            //         std::cout << node->name << ", ";
            // }
            // std::cout << std::endl;
        }

        // Run a single FM pass
        PartitionResult pass_result = run_fm_pass(); // run_fm_pass updates current partition state

        // std::cout << "[FM] Cut size after pass " << pass << " (best in pass): " << pass_result.cut_size << std::endl;

        // Check if the best result from this pass is better than the overall best
        if (pass_result.cut_size < best_result_overall.cut_size)
        { // Use the partition state *after* the pass rollback
            std::cout << "[FM] Pass " << pass << " found improvement. Overall best cut: " << pass_result.cut_size << std::endl;
            best_result_overall.cut_size = pass_result.cut_size;
            improvement_found = true;

            // Store the improved partition state
            best_partition_map.clear();
            for (Node *node : current_nodes)
            {
                best_partition_map[node] = node->partition_id;
            }
            // Also store current areas A and B if needed in result
            // best_result_overall.area_A = current_area_A;
            // best_result_overall.area_B = current_area_B;
        }
        else
        {
            // std::cout << "[FM] Pass " << pass << " did not improve overall best cut size (" << best_result_overall.cut_size << "). Stopping passes." << std::endl;
            // improvement_found remains false, loop will terminate
        }
    }

    // std::cout << "\n[FM] FM process completed. Best overall cut size: " << best_result_overall.cut_size << std::endl;

    // Restore the best partition found across all passes
    if (!best_partition_map.empty())
    {
        // std::cout << "[FM] Restoring best partition found." << std::endl;
        current_area_A = 0.0;
        current_area_B = 0.0;
        for (Node *node : current_nodes)
        {
            node->partition_id = best_partition_map[node];
            if (node->partition_id == 0)
            {
                current_area_A += node->area;
            }
            else
            {
                current_area_B += node->area;
            }
        }
    }
    else
    {
        // std::cout << "[FM] Warning: Best partition map was empty. Initial partition might be used." << std::endl;
        // Recalculate areas based on the final state just in case
        current_area_A = 0.0;
        current_area_B = 0.0;
        for (Node *node : current_nodes)
        {
            if (node->partition_id == 0)
            {
                current_area_A += node->area;
            }
            else
            {
                current_area_B += node->area;
            }
        }
    }

    // std::cout << "[FM] Final Partition Areas: A = " << current_area_A << ", B = " << current_area_B << std::endl;

    // Populate the final result structure (assuming it needs the node lists)
    best_result_overall.partition_A.clear();
    best_result_overall.partition_B.clear();
    for (Node *node : current_nodes)
    {
        if (node->partition_id == 0)
        {
            best_result_overall.partition_A.push_back(node);
        }
        else
        {
            best_result_overall.partition_B.push_back(node);
        }
    }

    return best_result_overall;
}

// --- Private Method Implementations ---

void FMPartitioner::initialize_partition(const std::vector<Node *> &nodes, double target_balance /* Unused */)
{
    // std::cout << "[FM] Initializing partition for " << nodes.size() << " nodes..." << std::endl;

    // Reset state
    reset_fm_state();
    // std::cout << "--- State after reset_fm_state ---" << std::endl;
    // for (auto cell : circuit_ref.cell_map)
    // { // Use the 'nodes' parameter directly
    //     if (!cell.first.empty())
    //     { // Basic null check
    //         std::cout << "  Node: " << cell.second.name
    //                   << ", Partition ID: " << cell.second.partition_id
    //                   << ", Gain: " << cell.second.fm_gain
    //                   << ", Locked: " << (cell.second.fm_locked ? "true" : "false")
    //                   << std::endl;
    //     }
    // }
    // std::cout << "----------------------------------" << std::endl;

    // Use a defined imbalance ratio instead of target_balance for initial assignment
    const double ALLOWED_AREA_IMBALANCE_RATIO = 0.1; // Example ratio (10%)

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1); // Randomly suggest partition 0 or 1

    for (Node *node : nodes)
    {
        int side = dis(gen); // 0 or 1

        if (side == 0) // Try partition A (id=0)
        {
            double prospective_area_A = current_area_A + node->area;
            double total_area = prospective_area_A + current_area_B;
            double diff = std::abs(prospective_area_A - current_area_B);
            // Handle potential division by zero or very small total area
            double ratio = (total_area <= 1e-9) ? 0.0 : (diff / total_area);

            if (ratio <= ALLOWED_AREA_IMBALANCE_RATIO)
            {
                node->partition_id = 0; // Assign to A
                current_area_A = prospective_area_A;
            }
            else
            {
                node->partition_id = 1; // Assign to B
                current_area_B += node->area;
            }
        }
        else // Try partition B (id=1)
        {
            double prospective_area_B = current_area_B + node->area;
            double total_area = current_area_A + prospective_area_B;
            double diff = std::abs(current_area_A - prospective_area_B);
            // Handle potential division by zero or very small total area
            double ratio = (total_area <= 1e-9) ? 0.0 : (diff / total_area);

            if (ratio <= ALLOWED_AREA_IMBALANCE_RATIO)
            {
                node->partition_id = 1; // Assign to B
                current_area_B = prospective_area_B;
            }
            else
            {
                node->partition_id = 0; // Assign to A
                current_area_A += node->area;
            }
        }
        // Note: node->fm_locked is reset in reset_fm_state()
    }

    // Print initial partition distribution (Optional)
    std::cout << "[FM] Initial Partition Areas: A = " << current_area_A << ", B = " << current_area_B << std::endl;
    long countA = 0, countB = 0;
    for (Node *node : nodes)
    {
        if (node->partition_id == 0)
            countA++;
        else
            countB++;
    }
    // std::cout << "[FM] Initial Partition Counts: A = " << countA << ", B = " << countB << std::endl;

    // std::cout << "--- State after random partition on currentNodes intialization ---" << std::endl;
    // for (auto cell : circuit_ref.cell_map)
    // { // Use the 'nodes' parameter directly
    //     if (!cell.first.empty())
    //     { // Basic null check
    //         std::cout << "  Node: " << cell.second.name
    //                   << ", Partition ID: " << cell.second.partition_id
    //                   << ", Gain: " << cell.second.fm_gain
    //                   << ", Locked: " << (cell.second.fm_locked ? "true" : "false")
    //                   << std::endl;
    //     }
    // }
    // std::cout << "----------------------------------" << std::endl;
}

void FMPartitioner::calculate_initial_gains()
{
    // std::cout << "[FM] Calculating initial gains..." << std::endl;

    // Reset gains for nodes to be partitioned
    for (Node *node_ptr : current_nodes)
    {
        node_ptr->fm_gain = 0;
    }

    // Reset net counts
    for (auto &net_pair : circuit_ref.net_map)
    {
        net_pair.second.partition_A_count = 0;
        net_pair.second.partition_B_count = 0;
    }

    // Calculate initial net distribution considering all nodes (including fixed ones)
    for (auto &node_pair : circuit_ref.cell_map) // Iterate with auto& for non-const access
    {
        Node &node = node_pair.second; // Get a non-const reference
        for (Net *net : node.nets)     // Access via reference
        {
            if (node.partition_id == 0) // Access via reference
            {
                net->partition_A_count++;
            }
            else if (node.partition_id == 1) // Assuming partition_id is 1 for B
            {
                net->partition_B_count++;
            }
        }
    }

    // Calculate gains based on net distribution for the nodes we are partitioning
    std::cout << "Gains for current Nodes" << std::endl;
    for (Node *node : current_nodes)
    {
        int gain = 0;
        int from_partition = node->partition_id;
        for (Net *net : node->nets)
        {
            int F_count = (from_partition == 0) ? net->partition_A_count : net->partition_B_count;
            int T_count = (from_partition == 0) ? net->partition_B_count : net->partition_A_count;

            if (F_count == 1)
                gain++;
            if (T_count == 0)
                gain--;
        }
        node->fm_gain = gain;
    }

    std::ofstream initial_gains_file("initial_gains.txt");
    for (Node *node : current_nodes)
    {
        initial_gains_file << node->name << "->" << node->fm_gain << ", partition_id->" << node->partition_id << std::endl;
    }
    initial_gains_file.close();
}

void FMPartitioner::build_gain_buckets()
{
    std::cout << "[FM] Building gain buckets..." << std::endl;

    // Clear existing buckets
    gain_bucket.clear();
    max_gain_index = -1; // Use -1 to indicate not yet found
    min_gain_index = std::numeric_limits<int>::max();

    // Find max potential gain (p_max)
    p_max = calculate_pmax(current_nodes);
    // Ensure bucket size is at least 1, even if p_max is 0
    std::cout << "P_MAX: " << p_max << std::endl;
    int bucket_size = std::max(1, 2 * p_max + 1);
    gain_bucket.resize(bucket_size);

    // Add nodes to appropriate buckets
    for (Node *node : current_nodes)
    {
        if (!node->fm_locked)
        {
            int bucket_index = node->fm_gain + p_max;

            // Clamp bucket_index to valid range [0, bucket_size - 1]
            bucket_index = std::max(0, std::min(bucket_index, bucket_size - 1));

            gain_bucket[bucket_index].push_front(node);
            node->bucket_iterator = gain_bucket[bucket_index].begin(); // Store iterator in the node

            // Update tracked max/min gain indices that have non-empty buckets
            if (bucket_index > max_gain_index)
            {
                max_gain_index = bucket_index;
            }
            // Use '<=' for min_gain_index initialization case
            if (bucket_index < min_gain_index)
            {
                min_gain_index = bucket_index;
            }
        }
    }
    std::cout << "Gain Buckets==>" << std::endl;
    for (int i = 0; i < gain_bucket.size(); i++)
    {
        std::cout << "i = " << i << ", nodes=";
        for (auto itr : gain_bucket[i])
        {
            std::cout << "(" << itr->name << ", " << itr->fm_gain << ") -->";
        }
        std::cout << std::endl;
    }

    // If no nodes were added, max_gain_index might still be -1
    if (max_gain_index == -1)
    {
        min_gain_index = 0; // Or handle case with no movable nodes
        max_gain_index = 0;
        return; // Nothing more to do
    }

    // Adjust min/max gain indices if buckets are empty at the extremes
    while (max_gain_index >= 0 && gain_bucket[max_gain_index].empty())
    {
        max_gain_index--;
    }
    // Ensure max_gain_index doesn't go below 0 if all buckets become empty (unlikely here)
    max_gain_index = std::max(0, max_gain_index);

    while (min_gain_index < gain_bucket.size() && gain_bucket[min_gain_index].empty())
    {
        min_gain_index++;
    }
    // Ensure min_gain_index doesn't exceed the actual bucket size or max_gain_index
    min_gain_index = std::min(min_gain_index, (int)gain_bucket.size() - 1);
    min_gain_index = std::min(min_gain_index, max_gain_index);
}

void FMPartitioner::reset_fm_state()
{
    current_area_A = 0.0;
    current_area_B = 0.0;
    current_partition_map.clear();
    lock_status.clear();

    // Iterate through ALL nodes in the circuit
    for (auto &pair : circuit_ref.cell_map)
    {
        Node &node = pair.second; // Get a reference to the Node object
        node.partition_id = -1;
        node.fm_gain = 0;
        node.fm_locked = true; // Set locked to true for all nodes as requested
    }
    // The loop below iterating over current_nodes is no longer needed for resetting these specific fields
    // as it's now handled by iterating over cell_map.
    // for (Node *node : current_nodes)
    // {
    //     node->partition_id = -1; // Ensure it's -1
    //     node->fm_gain = 0;
    //     node->fm_locked = true; // Set locked to true as requested
    //     // node->partition_id = -1; // Already done above
    // }
}

PartitionResult FMPartitioner::run_fm_pass()
{
    // std::cout << "[FM] Running FM pass..." << std::endl;

    PartitionResult best_pass_result;
    best_pass_result.cut_size = calculate_cut_size();
    best_pass_result.partition_A = {}; // TODO: Store initial partition state if needed for rollback
    best_pass_result.partition_B = {};

    int current_cut_size = best_pass_result.cut_size;
    double initial_area_A = current_area_A;
    double initial_area_B = current_area_B;

    // Store move sequence and corresponding gains/cuts for rollback
    std::vector<Node *> move_sequence;
    std::vector<int> gain_at_move;
    std::vector<int> cut_size_after_move;
    std::vector<double> area_A_after_move;
    std::vector<double> area_B_after_move;

    // Reset locks for the nodes we are partitioning
    for (Node *node : current_nodes)
    {
        node->fm_locked = false;
    }

    calculate_initial_gains();
    build_gain_buckets();

    int best_cut_in_pass = current_cut_size;
    int best_move_index = -1; // Index in move_sequence corresponding to best_cut_in_pass

    // Iterate while there are unlocked nodes with valid gains
    for (int move_count = 0; move_count < current_nodes.size(); ++move_count)
    {
        Node *node_to_move = select_best_node_to_move();

        if (!node_to_move)
        {
            // std::cout << "[FM] No feasible node found to move. Ending pass." << std::endl;
            break; // No more valid moves possible
        }

        int gain_of_move = node_to_move->fm_gain;

        // Move the node (updates gains of neighbors)
        move_node(node_to_move);

        // Record the move
        move_sequence.push_back(node_to_move);
        gain_at_move.push_back(gain_of_move);

        // Update cut size
        current_cut_size -= gain_of_move;
        cut_size_after_move.push_back(current_cut_size);
        area_A_after_move.push_back(current_area_A);
        area_B_after_move.push_back(current_area_B);

        // Track the best cut size seen in this pass
        if (current_cut_size < best_cut_in_pass)
        {
            best_cut_in_pass = current_cut_size;
            best_move_index = move_count;
        }
    }

    // std::cout << "[FM] Pass completed. Best cut in pass: " << best_cut_in_pass << " at move " << (best_move_index + 1) << std::endl;

    // Roll back moves to the state that gave the best cut size
    if (best_move_index != -1 && best_cut_in_pass < best_pass_result.cut_size)
    { // Only update if we found an improvement within the pass
        best_pass_result.cut_size = best_cut_in_pass;
        // Roll back from the end to the best state
        for (int i = move_sequence.size() - 1; i > best_move_index; --i)
        {
            Node *node_to_rollback = move_sequence[i];
            int original_partition = (node_to_rollback->partition_id == 0) ? 1 : 0;
            if (original_partition == 0)
            { // Moved to B, rollback to A
                current_area_B -= node_to_rollback->area;
                current_area_A += node_to_rollback->area;
            }
            else
            { // Moved to A, rollback to B
                current_area_A -= node_to_rollback->area;
                current_area_B += node_to_rollback->area;
            }
            node_to_rollback->partition_id = original_partition;
        }
        // Final area after rollback
        current_area_A = area_A_after_move[best_move_index];
        current_area_B = area_B_after_move[best_move_index];
    }
    else
    {
        // No improvement in this pass, revert all moves
        // std::cout << "[FM] No improvement found or best was initial state. Reverting all moves." << std::endl;
        for (int i = move_sequence.size() - 1; i >= 0; --i)
        {
            Node *node_to_rollback = move_sequence[i];
            int original_partition = (node_to_rollback->partition_id == 0) ? 1 : 0;
            if (original_partition == 0)
            { // Moved to B, rollback to A
                current_area_B -= node_to_rollback->area;
                current_area_A += node_to_rollback->area;
            }
            else
            { // Moved to A, rollback to B
                current_area_A -= node_to_rollback->area;
                current_area_B += node_to_rollback->area;
            }
            node_to_rollback->partition_id = original_partition;
        }
        // Reset area to initial state for this pass
        current_area_A = initial_area_A;
        current_area_B = initial_area_B;
    }

    // Populate PartitionResult with the final partition state after rollback
    // best_pass_result.partition_A.clear();
    // best_pass_result.partition_B.clear();
    // for(Node* node : current_nodes) {
    //     if (node->partition_id == 0) {
    //         best_pass_result.partition_A.push_back(node);
    //     } else {
    //         best_pass_result.partition_B.push_back(node);
    //     }
    // }

    return best_pass_result;
}

Node *FMPartitioner::select_best_node_to_move()
{
    // Iterate from the highest gain bucket downwards
    for (int current_index = max_gain_index; current_index >= min_gain_index; --current_index)
    {
        if (current_index < 0 || current_index >= gain_bucket.size())
            continue; // Boundary check

        auto &bucketList = gain_bucket[current_index];
        for (auto it = bucketList.begin(); it != bucketList.end(); /* manual increment */)
        {
            Node *candidate_node = *it;

            if (candidate_node->fm_locked)
            {
                ++it; // Skip locked nodes
                continue;
            }

            // Check if moving this node maintains area balance
            if (check_area_balance(candidate_node))
            {
                // Node is feasible, remove from bucket (will be locked in move_node)
                // Erasing invalidates iterator 'it', so get next before erasing if needed
                // but since we return immediately, just erase is fine.
                // bucketList.erase(it); // Remove here or after successful move?
                // Update max/min gain index if this bucket becomes empty
                // if (bucketList.empty()) {
                //    // Adjust max_gain_index or min_gain_index (complex logic, handle in update_single_node_gain)
                // }
                return candidate_node; // Found the best feasible node
            }
            else
            {
                ++it; // Check next node in the same bucket
            }
        }
    }

    // No feasible node found in any bucket
    return nullptr;
}

bool FMPartitioner::check_area_balance(Node *node_to_move)
{
    double prospective_area_A = current_area_A;
    double prospective_area_B = current_area_B;

    if (node_to_move->partition_id == 0)
    {
        // Moving A -> B
        prospective_area_A -= node_to_move->area;
        prospective_area_B += node_to_move->area;
    }
    else
    {
        // Moving B -> A
        prospective_area_B -= node_to_move->area;
        prospective_area_A += node_to_move->area;
    }

    double total_partitioned_area = current_area_A + current_area_B;
    // Avoid division by zero if total area is zero
    if (total_partitioned_area <= 0.0)
        return true;

    // Check if the prospective areas are within the allowed range
    bool balance_ok = (prospective_area_A >= current_target_area_A && prospective_area_A <= current_target_area_B) &&
                      (prospective_area_B >= current_target_area_A && prospective_area_B <= current_target_area_B);

    // Alternative check: Check total balance factor range
    // double balance_factor = prospective_area_A / total_partitioned_area;
    // bool balance_ok = (balance_factor >= balance_factor_min && balance_factor <= balance_factor_max);

    return balance_ok;
}

void FMPartitioner::move_node(Node *node_to_move)
{
    int from_partition = node_to_move->partition_id;
    int to_partition = (from_partition == 0) ? 1 : 0;

    // Update gains of neighbors *before* changing the partition counts
    update_gains_after_move(node_to_move);

    // Update areas
    if (from_partition == 0)
    { // Moving A -> B
        current_area_A -= node_to_move->area;
        current_area_B += node_to_move->area;
    }
    else
    { // Moving B -> A
        current_area_B -= node_to_move->area;
        current_area_A += node_to_move->area;
    }

    // Update net partition counts for nets connected to the moved node
    for (Net *net : node_to_move->nets)
    {
        if (from_partition == 0)
        {
            net->partition_A_count--;
            net->partition_B_count++;
        }
        else
        {
            net->partition_B_count--;
            net->partition_A_count++;
        }
    }

    // Change partition ID and lock the node
    node_to_move->partition_id = to_partition;
    node_to_move->fm_locked = true;

    // Use the stored iterator for O(1) removal
    if (node_to_move->bucket_iterator != std::list<Node *>::iterator())
    { // Check validity
        gain_bucket[node_to_move->fm_gain + p_max].erase(node_to_move->bucket_iterator);
    }
    else
    {
        // Fallback or warning
        // std::cerr << "[FM WARN] Node " << node_to_move->name << " had invalid bucket iterator during move_node. Trying list remove." << std::endl;
        gain_bucket[node_to_move->fm_gain + p_max].remove(node_to_move); // O(N) fallback
    }

    // We might need to update max/min gain indices if the bucket becomes empty
    // This logic is integrated into update_single_node_gain for efficiency
}

void FMPartitioner::update_gains_after_move(Node *moved_node)
{
    int from_partition = moved_node->partition_id;
    int to_partition = (from_partition == 0) ? 1 : 0;

    // Lock the moved node temporarily during updates to prevent self-update
    moved_node->fm_locked = true; // Already locked in move_node, but defensive

    for (Net *net : moved_node->nets)
    {
        // Cache counts before the hypothetical move
        int original_F_count = (from_partition == 0) ? net->partition_A_count : net->partition_B_count;
        int original_T_count = (from_partition == 0) ? net->partition_B_count : net->partition_A_count;

        // Check critical nets before the move
        if (original_T_count == 0)
        {
            // Moving the only node from T would make the net non-critical
            for (Node *neighbor_node : net->nodes)
            {
                if (!neighbor_node->fm_locked && neighbor_node != moved_node)
                {
                    update_single_node_gain(neighbor_node, +1);
                }
            }
        }
        else if (original_T_count == 1)
        {
            // Moving the node makes the only node in T critical
            for (Node *neighbor_node : net->nodes)
            {
                if (!neighbor_node->fm_locked && neighbor_node->partition_id == to_partition && neighbor_node != moved_node)
                {
                    update_single_node_gain(neighbor_node, -1);
                    // Optimization: Only one such neighbor exists in partition T
                    // break; // Careful if multiple nodes share the net
                }
            }
        }

        // --- Simulate the move by adjusting counts hypothetically ---
        // Decrement count in 'from' partition
        if (from_partition == 0)
            net->partition_A_count--;
        else
            net->partition_B_count--;
        // Increment count in 'to' partition
        if (to_partition == 0)
            net->partition_A_count++;
        else
            net->partition_B_count++;
        // --- Counts now reflect the state *after* the move ---

        // Cache counts after the hypothetical move
        int new_F_count = (from_partition == 0) ? net->partition_A_count : net->partition_B_count;
        int new_T_count = (from_partition == 0) ? net->partition_B_count : net->partition_A_count;

        // Check critical nets after the move
        if (new_F_count == 0)
        {
            // Moving this node made the net non-critical for the F side
            for (Node *neighbor_node : net->nodes)
            {
                if (!neighbor_node->fm_locked && neighbor_node != moved_node)
                {
                    update_single_node_gain(neighbor_node, -1);
                }
            }
        }
        else if (new_F_count == 1)
        {
            // Moving this node made the single node remaining in F critical
            for (Node *neighbor_node : net->nodes)
            {
                if (!neighbor_node->fm_locked && neighbor_node->partition_id == from_partition && neighbor_node != moved_node)
                {
                    update_single_node_gain(neighbor_node, +1);
                    // Optimization: Only one such neighbor exists in partition F
                    // break; // Careful if multiple nodes share the net
                }
            }
        }

        // --- Revert counts - actual update happens in move_node --- IS THIS NEEDED?
        // The update logic depends on counts *before* move. The actual counts are updated *after* this function in move_node.
        if (from_partition == 0)
            net->partition_A_count++;
        else
            net->partition_B_count++;
        if (to_partition == 0)
            net->partition_A_count--;
        else
            net->partition_B_count--;
        // --- Counts are back to state before the move --- Seems correct
    }

    // Unlock the node if necessary (it will be locked again in move_node)
    // moved_node->fm_locked = false; // Not needed, move_node locks it.
}

void FMPartitioner::update_single_node_gain(Node *node, int delta_gain)
{
    if (node->fm_locked)
        return; // Should not happen if called correctly

    int old_gain = node->fm_gain;
    int new_gain = old_gain + delta_gain;
    int old_bucket_index = old_gain + p_max;
    int new_bucket_index = new_gain + p_max;

    // Clamp indices to valid range
    old_bucket_index = std::max(0, std::min(old_bucket_index, (int)gain_bucket.size() - 1));
    new_bucket_index = std::max(0, std::min(new_bucket_index, (int)gain_bucket.size() - 1));

    // Remove node from the old bucket using its stored iterator
    if (node->bucket_iterator != std::list<Node *>::iterator())
    { // Check validity
        gain_bucket[old_bucket_index].erase(node->bucket_iterator);
    }
    else
    {
        // Fallback or warning
        // std::cerr << "[FM WARN] Node " << node->name << " had invalid bucket iterator during gain update (old bucket). Trying list remove." << std::endl;
        gain_bucket[old_bucket_index].remove(node); // O(N) fallback
    }

    // Add node to the new bucket
    gain_bucket[new_bucket_index].push_front(node);
    node->bucket_iterator = gain_bucket[new_bucket_index].begin(); // Store the new iterator in the node

    // Update node's gain value
    node->fm_gain = new_gain;

    // Update max/min gain indices if necessary
    // If the old bucket became empty and it was an extreme
    if (gain_bucket[old_bucket_index].empty())
    {
        if (old_bucket_index == max_gain_index)
        {
            while (max_gain_index >= 0 && gain_bucket[max_gain_index].empty())
            {
                max_gain_index--;
            }
            max_gain_index = std::max(0, max_gain_index); // Don't go below 0
        }
        if (old_bucket_index == min_gain_index)
        {
            while (min_gain_index < gain_bucket.size() && gain_bucket[min_gain_index].empty())
            {
                min_gain_index++;
            }
            min_gain_index = std::min((int)gain_bucket.size() - 1, min_gain_index); // Don't exceed size
            min_gain_index = std::min(min_gain_index, max_gain_index);              // Ensure min <= max
        }
    }

    // If the new bucket index extends the range
    if (new_bucket_index > max_gain_index)
    {
        max_gain_index = new_bucket_index;
    }
    if (new_bucket_index < min_gain_index)
    {
        min_gain_index = new_bucket_index;
    }
}

int FMPartitioner::calculate_cut_size()
{
    int cut = 0;
    // Iterate through all nets in the circuit
    for (const auto &net_pair : circuit_ref.net_map)
    {
        const Net &net = net_pair.second;
        bool has_A = false;
        bool has_B = false;
        // Check the partition ID of each node connected to the net
        for (const Node *node : net.nodes)
        {
            if (node->partition_id == -1)
                continue;
            if (node->partition_id == 0)
            {
                has_A = true;
            }
            else if (node->partition_id == 1)
            {
                has_B = true;
            }
            // Optimization: if both partitions are found, no need to check further
            if (has_A && has_B)
            {
                break;
            }
        }
        // If nodes from both partitions are present, the net is cut
        if (has_A && has_B)
        {
            cut++;
        }
    }
    return cut;
}

double FMPartitioner::calculate_total_area(const std::vector<Node *> &nodes)
{
    // srand(time(0)); // Should not be here - call once in main or constructor if needed
    double total = 0.0;
    for (Node *node : nodes)
    {
        // Only sum area for nodes we are partitioning
        total += node->area;
    }
    return total;
}

int FMPartitioner::calculate_pmax(const std::vector<Node *> &nodes)
{
    int max_degree = 0;
    for (Node *node : nodes)
    {
        // Consider only nodes being partitioned for p_max calculation relevant to gains
        if (node)
        {
            max_degree = std::max(max_degree, static_cast<int>(node->nets.size()));
        }
    }
    // p_max should be at least 1 to avoid issues with bucket indexing if all nodes have degree 0
    return std::max(1, max_degree);
}