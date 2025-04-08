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

    balance_factor_min_d = balance_factor_min;
    balance_factor_max_d = balance_factor_max;

    // std::cout << "[FM] Target Area Range: [" << current_target_area_A << ", " << current_target_area_B << "]" << std::endl;

    // Initialize the partition randomly for the nodes to be partitioned
    initialize_partition(current_nodes, (balance_factor_min + balance_factor_max) / 2.0); // Use average balance for initial random split

    // Store the best overall partition found
    PartitionResult best_result_overall;
    best_result_overall.cut_size = std::numeric_limits<int>::max();
    std::unordered_map<Node *, int> best_partition_map;

    int max_passes = 2; // Maximum number of passes
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
            // std::cout << "[FM] Pass " << pass << " found improvement. Overall best cut: " << pass_result.cut_size << std::endl;
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
    // std::cout << "[FM] Initial Partition Areas: A = " << current_area_A << ", B = " << current_area_B << std::endl;
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

void FMPartitioner::initialize_gains_and_buckets()
{
    // std::cout << "[FM] Initializing gains and building buckets..." << std::endl;

    // 1. Reset state
    gain_bucket.clear();
    p_max = 0;
    max_gain_index = std::numeric_limits<int>::min();
    min_gain_index = std::numeric_limits<int>::max();

    // Reset gains for nodes to be partitioned (unlocking happens in initialize_partition)
    for (Node *node_ptr : current_nodes)
    {
        node_ptr->fm_gain = 0;
        node_ptr->fm_locked = false;
    }

    // 2. Reset Net Counts
    for (auto &net_pair : circuit_ref.net_map)
    {
        net_pair.second.partition_A_count = 0;
        net_pair.second.partition_B_count = 0;
    }

    // 3. Calculate Initial Net Distribution based on the current partition state
    for (auto const &[node_name, node] : circuit_ref.cell_map)
    {
        p_max = std::max(p_max, static_cast<int>(node.nets.size()));
        for (Net *net : node.nets)
        {
            if (node.partition_id == 0)
                net->partition_A_count++;
            else if (node.partition_id == 1)
                net->partition_B_count++;
        }
    }

    // 5. Resize Gain Buckets
    int bucket_size = 2 * p_max + 1;
    gain_bucket.resize(bucket_size);

    // 6. Calculate Gains and Populate Buckets for current_nodes
    bool nodes_added = false;
    max_gain_index = 0;               // Reset effective max index
    min_gain_index = bucket_size - 1; // Reset effective min index

    for (Node *node : current_nodes)
    {
        int gain = 0;
        for (Net *net : node->nets)
        {
            int F_count = (node->partition_id == 0) ? net->partition_A_count : net->partition_B_count;
            int T_count = (node->partition_id == 0) ? net->partition_B_count : net->partition_A_count;

            if (F_count == 1)
                gain++;
            if (T_count == 0)
                gain--;
        }
        node->fm_gain = gain;
        int bucket_index = gain + p_max;
        gain_bucket[bucket_index].push_front(node);
        node->bucket_iterator = gain_bucket[bucket_index].begin();
        max_gain_index = std::max(max_gain_index, bucket_index);
        min_gain_index = std::min(min_gain_index, bucket_index);
    }

    // printGainBucket();
}

void FMPartitioner::printGainBucket()
{
    std::cout << "[INFO] Gain Buckets " << std::endl;
    for (int i = 0; i < gain_bucket.size(); i++)
    {
        std::cout << "Gain Bucket[" << i << "](" << i - p_max << "): ";
        for (auto &cellId : gain_bucket[i])
        {
            std::cout << cellId->name << " ";
        }
        std::cout << std::endl;
    }
}

void FMPartitioner::reset_fm_state()
{
    current_area_A = 0.0;
    current_area_B = 0.0;

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
    int THRESHOLD = 470000;
    PartitionResult best_pass_result;
    best_pass_result.cut_size = calculate_cut_size();
    best_pass_result.partition_A = {}; // TODO: Store initial partition state if needed for rollback
    best_pass_result.partition_B = {};

    int current_cut_size = best_pass_result.cut_size;
    double initial_area_A = current_area_A;
    double initial_area_B = current_area_B;

    // Store move sequence and corresponding gains/cuts for rollback
    std::vector<Node *> move_sequence;
    std::vector<int> cut_sizes;
    std::vector<double> area_A_history;
    std::vector<double> area_B_history;

    // Initialize tracking vectors with initial state
    cut_sizes.push_back(current_cut_size);
    area_A_history.push_back(current_area_A);
    area_B_history.push_back(current_area_B);

    initialize_gains_and_buckets();

    int best_cut_in_pass = current_cut_size;
    int best_move_index = -1; // Index in move_sequence corresponding to best_cut_in_pass

    int iter = 0;

    // Iterate while there are unlocked nodes with valid gains
    while (max_gain_index >= min_gain_index)
    {
        iter++;

        int currMaxGain = max_gain_index - p_max;
        auto &bucketList = gain_bucket[max_gain_index];

        if (bucketList.empty())
        {
            max_gain_index--;
            continue;
        }

        auto &targetNode = bucketList.back();
        bucketList.pop_back();
        if (bucketList.empty())
            max_gain_index--;

        if (targetNode->fm_locked)
            continue;

        if (!check_area_balance(targetNode))
        {
            targetNode->fm_locked = true;
            continue;
        }

        // Store the node being moved
        move_sequence.push_back(targetNode);

        if (iter > THRESHOLD)
        {
            std::cout << "Start: " << std::endl;
        }

        // Update gains and move the node
        update_gains_after_move(targetNode);
        if (iter > THRESHOLD)
        {
            std::cout << "End: " << std::endl;
        }

        targetNode->fm_locked = true;
        int oldPartition = targetNode->partition_id;
        targetNode->partition_id = oldPartition == 0 ? 1 : 0;

        // Update current cut size and areas
        // current_cut_size = calculate_cut_size();
        current_cut_size = current_cut_size - targetNode->fm_gain;
        cut_sizes.push_back(current_cut_size);
        area_A_history.push_back(current_area_A);
        area_B_history.push_back(current_area_B);

        // Check if this is the best cut so far
        if (current_cut_size < best_cut_in_pass)
        {
            best_cut_in_pass = current_cut_size;
            best_move_index = move_sequence.size() - 1;
        }

        if (current_cut_size > 100000 && current_cut_size > 2 * best_cut_in_pass)
        {
            break;
        }

        if (iter > THRESHOLD)
        {
            std::cout << "End2: " << std::endl;
        }

        if (iter > THRESHOLD)
        {
            // if (iter % 100 == 0)
            {

                std::cout << "Gainbucket: " << gain_bucket.size() << std::endl;
                std::cout << "current: " << current_cut_size << std::endl;
                std::cout << "BucketList: " << bucketList.size() << std::endl;
            }
        }
    }

    // If we found a better cut, rollback to that state
    if (best_move_index >= 0)
    {
        // std::cout << "Cutsizes: " << cut_sizes.size() << " best_move_index: " << best_move_index << std::endl;
        // for (auto cut_size : cut_sizes)
        // {
        //     std::cout << cut_size << " , ";
        // }
        // std::cout << std::endl;
        // std::cout << "Move Sequence(" << move_sequence.size() << "): ";
        // for (auto node : move_sequence)
        // {
        //     std::cout << node->name << " , ";
        // }
        // std::cout << std::endl;
        rollback_to_best_cut(move_sequence, cut_sizes, area_A_history, area_B_history, best_move_index, best_pass_result);
        best_pass_result.cut_size = best_cut_in_pass;
    }

    // std::cout << "[FM] Pass completed. Best cut in pass: " << best_cut_in_pass
    //           << " at move&cut " << (best_move_index + 1) << "&" << cut_sizes[best_move_index + 1] << std::endl;

    return best_pass_result;
}

bool FMPartitioner::check_area_balance(Node *node_to_move)
{
    double prospective_area_A = current_area_A;
    double prospective_area_B = current_area_B;

    if (node_to_move->partition_id == 0)
    {
        prospective_area_A -= node_to_move->area;
        prospective_area_B += node_to_move->area;
    }
    else
    {
        prospective_area_B -= node_to_move->area;
        prospective_area_A += node_to_move->area;
    }

    double total_area = current_area_A + current_area_B;
    if (total_area <= 0.0)
        return true;

    double balance_factor = prospective_area_A / total_area;
    return (balance_factor >= balance_factor_min_d && balance_factor <= balance_factor_max_d);
}

// void FMPartitioner::move_node(Node *node_to_move)
// {
//     int from_partition = node_to_move->partition_id;
//     int to_partition = (from_partition == 0) ? 1 : 0;

//     // Update gains of neighbors *before* changing the partition counts
//     update_gains_after_move(node_to_move);

//     // Update areas
//     if (from_partition == 0)
//     { // Moving A -> B
//         current_area_A -= node_to_move->area;
//         current_area_B += node_to_move->area;
//     }
//     else
//     { // Moving B -> A
//         current_area_B -= node_to_move->area;
//         current_area_A += node_to_move->area;
//     }

//     // Update net partition counts for nets connected to the moved node
//     for (Net *net : node_to_move->nets)
//     {
//         if (from_partition == 0)
//         {
//             net->partition_A_count--;
//             net->partition_B_count++;
//         }
//         else
//         {
//             net->partition_B_count--;
//             net->partition_A_count++;
//         }
//     }

//     // Change partition ID and lock the node
//     node_to_move->partition_id = to_partition;
//     node_to_move->fm_locked = true;

//     // Use the stored iterator for O(1) removal
//     if (node_to_move->bucket_iterator != std::list<Node *>::iterator())
//     { // Check validity
//         gain_bucket[node_to_move->fm_gain + p_max].erase(node_to_move->bucket_iterator);
//     }
//     else
//     {
//         // Fallback or warning
//         // std::cerr << "[FM WARN] Node " << node_to_move->name << " had invalid bucket iterator during move_node. Trying list remove." << std::endl;
//         gain_bucket[node_to_move->fm_gain + p_max].remove(node_to_move); // O(N) fallback
//     }

//     // We might need to update max/min gain indices if the bucket becomes empty
//     // This logic is integrated into update_single_node_gain for efficiency
// }

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
    moved_node->fm_locked = true; // Not needed, move_node locks it.
}

void FMPartitioner::update_single_node_gain(Node *node, int delta_gain)
{
    if (node->fm_locked)
        return;

    int old_gain = node->fm_gain;
    int new_gain = old_gain + delta_gain;
    int old_bucket_index = old_gain + p_max;
    int new_bucket_index = new_gain + p_max;

    // Remove node from the old bucket
    gain_bucket[old_bucket_index].erase(node->bucket_iterator);

    // Add node to the new bucket
    gain_bucket[new_bucket_index].push_front(node);
    node->bucket_iterator = gain_bucket[new_bucket_index].begin();

    // Update node's gain
    node->fm_gain = new_gain;

    // Update max/min gain indices only if necessary
    if (gain_bucket[old_bucket_index].empty())
    {
        if (old_bucket_index == max_gain_index)
        {
            while (max_gain_index >= 0 && gain_bucket[max_gain_index].empty())
                max_gain_index--;
        }
        if (old_bucket_index == min_gain_index)
        {
            while (min_gain_index < gain_bucket.size() && gain_bucket[min_gain_index].empty())
                min_gain_index++;
        }
    }
    max_gain_index = std::max(max_gain_index, new_bucket_index);
    min_gain_index = std::min(min_gain_index, new_bucket_index);
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
    double total = 0.0;
    for (Node *node : nodes)
    {
        total += node->area;
    }
    return total;
}

void FMPartitioner::rollback_to_best_cut(const std::vector<Node *> &move_sequence,
                                         const std::vector<int> &cut_sizes,
                                         const std::vector<double> &area_A_history,
                                         const std::vector<double> &area_B_history,
                                         int best_move_index, PartitionResult best_pass_result)
{
    // Validate inputs and show which validations fail
    bool validation_failed = false;
    std::string error_msg = "[FM] Invalid rollback parameters:";

    if (move_sequence.empty())
    {
        validation_failed = true;
        error_msg += "\n  - Move sequence is empty";
    }
    if (best_move_index < 0)
    {
        validation_failed = true;
        error_msg += "\n  - Best move index is negative";
    }
    if (best_move_index >= static_cast<int>(move_sequence.size()))
    {
        validation_failed = true;
        error_msg += "\n  - Best move index exceeds move sequence size";
    }
    if (cut_sizes.size() != move_sequence.size() + 1)
    {
        validation_failed = true;
        error_msg += "\n  - Cut sizes size (" + std::to_string(cut_sizes.size()) +
                     ") does not match expected size (" + std::to_string(move_sequence.size() + 1) + ")";
    }
    if (area_A_history.size() != move_sequence.size() + 1)
    {
        validation_failed = true;
        error_msg += "\n  - Area A history size (" + std::to_string(area_A_history.size()) +
                     ") does not match expected size (" + std::to_string(move_sequence.size() + 1) + ")";
    }
    if (area_B_history.size() != move_sequence.size() + 1)
    {
        validation_failed = true;
        error_msg += "\n  - Area B history size (" + std::to_string(area_B_history.size()) +
                     ") does not match expected size (" + std::to_string(move_sequence.size() + 1) + ")";
    }

    if (validation_failed)
    {
        std::cerr << error_msg << std::endl;
        return;
    }

    // Reset all nodes to their initial state (before any moves)
    for (Node *node : best_pass_result.partition_A)
    {
        node->partition_id = 0; // Reset to initial partition
        node->fm_locked = false;
    }

    // Reset all nodes to their initial state (before any moves)
    for (Node *node : best_pass_result.partition_B)
    {
        node->partition_id = 1; // Reset to initial partition
        node->fm_locked = false;
    }

    // Replay moves up to the best cut point
    for (int i = 0; i <= best_move_index; ++i)
    {
        Node *node = move_sequence[i];
        // Move the node to the opposite partition
        node->partition_id = (node->partition_id == 0) ? 1 : 0;
    }

    // Update the current areas to match the best cut state
    current_area_A = area_A_history[best_move_index + 1];
    current_area_B = area_B_history[best_move_index + 1];

    // Reinitialize gains and buckets for the best partition state
    initialize_gains_and_buckets();

    // std::cout << "[FM] Rolled back to best cut state with cut size: "
    //           << cut_sizes[best_move_index + 1] << "on move " << best_move_index + 1 << std::endl;
}