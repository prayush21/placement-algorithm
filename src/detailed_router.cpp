#include "detailed_router.h"
#include <iostream>  // For std::cout, std::cerr (debugging)
#include <algorithm> // For std::min, std::max, std::sort
#include <limits>    // For std::numeric_limits

namespace router
{

    DetailedRouter::DetailedRouter()
        : channel_length_(0), channel_capacity_(0), num_tracks_used_(0)
    {
        // Constructor: Initialize members
    }

    DetailedRouter::~DetailedRouter()
    {
        // Destructor: Clean up if necessary (e.g., if DetailedRouter allocated memory for NetChannelRepresentation objects)
        // However, based on current design, Channel owns NetChannelRepresentation objects if they are heap-allocated,
        // or they are part of Channel.nets directly if ChannelParser manages their lifecycle.
    }

    void DetailedRouter::initialize_net_representations(const std::vector<NetChannelRepresentation *> &nets)
    {
        std::cout << "Initializing net representations..." << std::endl;
        for (NetChannelRepresentation *net_repr : nets)
        {
            if (!net_repr || net_repr->pins.empty())
            {
                // std::cerr << "Warning: Net " << (net_repr ? net_repr->net_id_str : "UNKNOWN") << " has no pins or is null." << std::endl;
                if (net_repr)
                {
                    net_repr->left_most_x = std::numeric_limits<int>::max();
                    net_repr->right_most_x = std::numeric_limits<int>::min();
                }
                continue;
            }

            int min_x = std::numeric_limits<int>::max();
            int max_x = std::numeric_limits<int>::min();

            for (const auto &pin : net_repr->pins)
            {
                min_x = std::min(min_x, pin.x_coord);
                max_x = std::max(max_x, pin.x_coord);
            }
            net_repr->left_most_x = min_x;
            net_repr->right_most_x = max_x;
            net_repr->is_routed = false;
            net_repr->track_assignment = -1;
            net_repr->is_primary_segment = true;
            net_repr->current_segment_start_x = min_x;
            net_repr->current_segment_end_x = max_x;
        }
    }

    void DetailedRouter::build_vertical_constraint_graph(const std::vector<NetChannelRepresentation *> &routable_nets)
    {
        std::cout << "Building VCG..." << std::endl;
        vcg_above_.clear();
        vcg_below_.clear();

        for (int x = 0; x < channel_length_; ++x)
        {
            NetChannelRepresentation *top_pin_net = nullptr;
            NetChannelRepresentation *bottom_pin_net = nullptr;

            for (NetChannelRepresentation *net_repr : routable_nets)
            {
                for (const auto &pin : net_repr->pins)
                {
                    if (pin.x_coord == x)
                    {
                        if (pin.side == PinSide::TOP)
                        {
                            top_pin_net = net_repr;
                        }
                        else if (pin.side == PinSide::BOTTOM)
                        {
                            bottom_pin_net = net_repr;
                        }
                    }
                }
            }

            if (top_pin_net && bottom_pin_net && top_pin_net != bottom_pin_net)
            {
                vcg_above_[bottom_pin_net].insert(top_pin_net);
                vcg_below_[top_pin_net].insert(bottom_pin_net);
            }
        }
    }

    // void DetailedRouter::calculate_vcg_in_degrees() // Remove this definition
    // {
    //     std::cout << "Calculating VCG in-degrees..." << std::endl;
    //     vcg_in_degree_.clear();
    //     // ... rest of the commented out function ...
    // }

    NetChannelRepresentation *DetailedRouter::select_next_net_for_lea(std::vector<NetChannelRepresentation *> &active_nets)
    {
        NetChannelRepresentation *best_net = nullptr;
        // Among candidates (active_nets whose left edge is <= current column and not yet routed):
        // Select one with VCG in-degree of 0 (from unrouted predecessors).
        // Tie-breaking: e.g., leftmost end point, shortest net, original net ID.
        int min_left_x = std::numeric_limits<int>::max();

        for (NetChannelRepresentation *net : active_nets)
        {
            if (net->is_routed)
                continue;

            bool has_unrouted_predecessors = false;
            if (vcg_below_.count(net))
            {
                for (const auto *predecessor_net : vcg_below_.at(net))
                {
                    if (!predecessor_net->is_routed)
                    {
                        has_unrouted_predecessors = true;
                        break;
                    }
                }
            }

            if (!has_unrouted_predecessors)
            {
                if (net->current_segment_start_x < min_left_x)
                {
                    min_left_x = net->current_segment_start_x;
                    best_net = net;
                }
                else if (net->current_segment_start_x == min_left_x)
                {
                    if (best_net == nullptr || net->current_segment_end_x < best_net->current_segment_end_x)
                    {
                        best_net = net;
                    }
                }
            }
        }
        return best_net;
    }

    bool DetailedRouter::assign_net_to_track(NetChannelRepresentation *net_repr, int track_idx, int segment_start_x, int segment_end_x)
    {
        if (track_idx < 0 || track_idx >= tracks_layout_.size())
        {
            // std::cerr << "Error: Invalid track index " << track_idx << std::endl;
            return false;
        }
        if (!tracks_layout_[track_idx].can_place(net_repr, segment_start_x, segment_end_x))
        {
            // std::cout << "Cannot place net " << net_repr->net_id_str << " on track " << track_idx << " due to overlap." << std::endl;
            return false;
        }

        tracks_layout_[track_idx].add_segment(net_repr, segment_start_x, segment_end_x);
        net_repr->track_assignment = track_idx;
        // For simple LEA without doglegs, is_routed would be set to true.
        // With doglegs, only the current segment might be considered 'done', and is_routed might mean fully done.
        // net_repr->is_routed = true; // This needs refinement for doglegs.
        // For now, assume one segment per net.
        if (net_repr->is_primary_segment && segment_end_x >= net_repr->right_most_x)
        {
            net_repr->is_routed = true; // Mark as fully routed if this segment covers the whole span
        }

        // std::cout << "  Assigned Net " << net_repr->net_id_str << " segment [" << segment_start_x << ", " << segment_end_x
        //           << "] to track " << track_idx << std::endl;
        return true;
    }

    bool DetailedRouter::left_edge_algorithm_route(const std::vector<NetChannelRepresentation *> &initial_routable_nets, Channel &channel_data_ref_for_length_capacity)
    {
        std::cout << "Applying Left-Edge Algorithm for channel " << channel_data_ref_for_length_capacity.id << "..." << std::endl;
        num_tracks_used_ = 0;
        tracks_layout_.clear();

        int current_max_tracks = channel_capacity_ > 0 ? channel_capacity_ : 1;
        tracks_layout_.resize(current_max_tracks, Track(0));
        for (int i = 0; i < current_max_tracks; ++i)
            tracks_layout_[i].track_idx = i;

        // Make a mutable copy of the initial_routable_nets to manage segments (original, merged, doglegs)
        std::vector<NetChannelRepresentation *> work_list_nets = initial_routable_nets;
        // Sort this work_list once, as new dogleg segments will be added and might need re-evaluation of order or specific handling.
        std::sort(work_list_nets.begin(), work_list_nets.end(), [](const NetChannelRepresentation *a, const NetChannelRepresentation *b)
                  {
                      if (a->current_segment_start_x != b->current_segment_start_x)
                      {
                          return a->current_segment_start_x < b->current_segment_start_x;
                      }
                      return a->net_id_str < b->net_id_str; // Tie-breaking
                  });

        std::vector<NetChannelRepresentation *> active_nets_list;            // Nets currently spanning column x
        std::vector<NetChannelRepresentation *> pending_new_dogleg_segments; // Doglegs created in current column sweep

        int net_idx_to_process = 0; // Index for iterating through work_list_nets

        for (int x = 0; x < channel_length_; ++x)
        {
            // Add nets starting at column x from the work_list to active_nets_list
            while (net_idx_to_process < work_list_nets.size() &&
                   work_list_nets[net_idx_to_process]->current_segment_start_x == x)
            {
                if (!work_list_nets[net_idx_to_process]->is_routed)
                { // Only add if not already fully routed
                    active_nets_list.push_back(work_list_nets[net_idx_to_process]);
                }
                net_idx_to_process++;
            }

            // Remove nets from active_nets_list if their current segment ends before column x or they got routed.
            active_nets_list.erase(
                std::remove_if(active_nets_list.begin(), active_nets_list.end(),
                               [x](NetChannelRepresentation *n)
                               {
                                   return n->is_routed || n->current_segment_end_x < x;
                               }),
                active_nets_list.end());

            // Add any newly created dogleg segments from the PREVIOUS column sweep if they start at x
            if (!pending_new_dogleg_segments.empty())
            {
                std::vector<NetChannelRepresentation *> still_pending;
                for (NetChannelRepresentation *seg : pending_new_dogleg_segments)
                {
                    if (seg->current_segment_start_x == x)
                    {
                        active_nets_list.push_back(seg);
                        // Add to main work_list as well for future global sort/processing if LEA iterates multiple times
                        // For single pass LEA, active_list is enough if sorted properly.
                        // Check if it's already in work_list_nets to avoid duplicates (if added from multiple places)
                        bool already_in_work_list = false;
                        for (const auto *wn : work_list_nets)
                            if (wn == seg)
                                already_in_work_list = true;
                        if (!already_in_work_list)
                            work_list_nets.push_back(seg);
                    }
                    else if (seg->current_segment_start_x > x)
                    {
                        still_pending.push_back(seg); // Keep for future columns
                    }
                    // Segments starting before x should have been processed or are errors
                }
                pending_new_dogleg_segments = still_pending;
                // Re-sort active_nets_list if new segments were added
                std::sort(active_nets_list.begin(), active_nets_list.end(), [](const NetChannelRepresentation *a, const NetChannelRepresentation *b)
                          {
                    if (a->current_segment_start_x != b->current_segment_start_x) return a->current_segment_start_x < b->current_segment_start_x;
                    if (a->current_segment_end_x != b->current_segment_end_x) return a->current_segment_end_x < b->current_segment_end_x;
                    return a->net_id_str < b->net_id_str; });
            }

            NetChannelRepresentation *net_to_route = nullptr;
            while ((net_to_route = select_next_net_for_lea(active_nets_list)) != nullptr)
            {
                bool placed_successfully = false;
                for (int t_idx = 0; t_idx < tracks_layout_.size(); ++t_idx)
                {
                    if (assign_net_to_track(net_to_route, t_idx, net_to_route->current_segment_start_x, net_to_route->current_segment_end_x))
                    {
                        placed_successfully = true;
                        num_tracks_used_ = std::max(num_tracks_used_, t_idx + 1);
                        NetChannelRepresentation *ultimate_parent = net_to_route->parent_representation ? net_to_route->parent_representation : net_to_route;
                        if (net_to_route->current_segment_end_x >= ultimate_parent->right_most_x)
                        {
                            ultimate_parent->is_routed = true;
                        }
                        net_to_route->is_routed = true;
                        break;
                    }
                }

                if (!placed_successfully)
                {
                    if (tracks_layout_.size() < channel_capacity_ || channel_capacity_ == 0)
                    {
                        tracks_layout_.emplace_back(tracks_layout_.size());
                        if (assign_net_to_track(net_to_route, tracks_layout_.size() - 1, net_to_route->current_segment_start_x, net_to_route->current_segment_end_x))
                        {
                            placed_successfully = true;
                            num_tracks_used_ = std::max(num_tracks_used_, (int)tracks_layout_.size());
                            NetChannelRepresentation *ultimate_parent = net_to_route->parent_representation ? net_to_route->parent_representation : net_to_route;
                            if (net_to_route->current_segment_end_x >= ultimate_parent->right_most_x)
                            {
                                ultimate_parent->is_routed = true;
                            }
                            net_to_route->is_routed = true;
                        }
                    }
                }

                if (!placed_successfully)
                {
                    std::cout << "  Net " << net_to_route->net_id_str << " segment [" << net_to_route->current_segment_start_x << "," << net_to_route->current_segment_end_x << "] could not be placed directly. Attempting dogleg..." << std::endl;
                    bool dogleg_created_this_attempt = false;
                    // Try to split at current column x (jog occurs at x, first part ends at x-1 or x, second starts at x or x+1)
                    // A better heuristic would be to try splitting at each pin of net_to_route, or just before a conflict.
                    // For now, we try to split such that the first part ends at `x` (the current column).
                    // The `attempt_dogleg_split` expects `proposed_dogleg_x` to be the column where the jog happens (end of first part).
                    int proposed_split_col = x;
                    if (proposed_split_col > net_to_route->current_segment_start_x && proposed_split_col < net_to_route->current_segment_end_x)
                    {
                        for (int t_idx = 0; t_idx < tracks_layout_.size(); ++t_idx)
                        {
                            NetChannelRepresentation *new_segment = attempt_dogleg_split(net_to_route, proposed_split_col, t_idx);
                            if (new_segment)
                            {
                                pending_new_dogleg_segments.push_back(new_segment);
                                num_tracks_used_ = std::max(num_tracks_used_, t_idx + 1);
                                // net_to_route (the first part) is now shorter. Its is_routed status is handled by attempt_dogleg_split/assign.
                                // If the first part of the dogleg is considered routed (e.g. it ends at its original right_most_x, unlikely for a useful dogleg)
                                // or if its current_segment_end_x makes it complete for *its new shortened span*.
                                // Check if the *first part* (net_to_route) is now considered routed for its (new) span.
                                // This is implicitly handled if net_to_route->current_segment_end_x (updated by dogleg_split) is >= its (possibly also updated) right_most_x.
                                // For now, rely on assign_net_to_track logic or dogleg_split to set is_routed for the first part if applicable.
                                // The crucial part is that the first segment (net_to_route) is now on a track.
                                net_to_route->is_routed = true; // The *first part* of the dogleg is now routed.

                                placed_successfully = true; // Overall placement for *this iteration* of net_to_route succeeded (by splitting).
                                dogleg_created_this_attempt = true;
                                std::cout << "    Dogleg successful for " << net_to_route->net_id_str << ". First part on track " << t_idx << ". New segment " << new_segment->net_id_str << " pending." << std::endl;
                                break;
                            }
                        }
                    }
                    else
                    {
                        std::cout << "    Cannot dogleg net " << net_to_route->net_id_str << " at column " << x << " (split point " << proposed_split_col << " invalid for segment ["
                                  << net_to_route->current_segment_start_x << "," << net_to_route->current_segment_end_x << "])." << std::endl;
                    }

                    if (!dogleg_created_this_attempt)
                    {
                        std::cerr << "    Dogleg attempt failed for net " << net_to_route->net_id_str << " segment." << std::endl;
                        net_to_route->is_routed = true;      // Mark as processed to avoid infinite loop
                        net_to_route->track_assignment = -2; // Indicate routing failure for this segment/net
                    }
                }

                if (placed_successfully || net_to_route->is_routed)
                {
                    active_nets_list.erase(std::remove(active_nets_list.begin(), active_nets_list.end(), net_to_route), active_nets_list.end());
                }
                else
                {
                    // If it wasn't placed, and wasn't marked routed (e.g. due to dogleg failure path setting is_routed = true)
                    // it means it's still in active_nets_list and select_next_net_for_lea might pick it again.
                    // This path should ideally not be hit if all failure modes set is_routed=true for the current segment.
                    // To be safe, remove it from active_nets_list if it could not be processed at all in this iteration.
                    std::cerr << "Warning: Net " << net_to_route->net_id_str << " was not processed and not marked routed/failed. Removing from active list to prevent loop." << std::endl;
                    active_nets_list.erase(std::remove(active_nets_list.begin(), active_nets_list.end(), net_to_route), active_nets_list.end());
                }
            }
        }

        // After all columns, add any remaining pending doglegs to the main work list if they weren't processed
        // This is more for multi-pass LEA. In single pass, they should have been added to active_nets_list when their start_x was reached.
        for (auto *seg : pending_new_dogleg_segments)
        {
            bool exists = false;
            for (auto *wl_net : work_list_nets)
                if (wl_net == seg)
                    exists = true;
            if (!exists)
                work_list_nets.push_back(seg);
        }
        pending_new_dogleg_segments.clear();

        // Final check for unrouted nets (original or merged, considering all their segments)
        bool all_original_nets_routed = true;
        for (NetChannelRepresentation *initial_net_repr : initial_routable_nets)
        {
            // An initial net (original or merged) is considered routed if its `is_routed` flag is true.
            // This flag should be set by the LEA logic when all its segments (if any doglegs) are placed
            // and cover the full span of this initial_net_repr.
            // The current logic sets `ultimate_parent->is_routed = true` when a segment ending at `ultimate_parent->right_most_x` is placed.
            if (!initial_net_repr->is_routed && !(initial_net_repr->pins.empty()))
            {
                // Check if any of its generated segments failed
                bool failed_segment_exists = false;
                for (const auto &gen_net_ptr : owned_generated_nets_)
                {
                    if (gen_net_ptr->parent_representation == initial_net_repr && gen_net_ptr->track_assignment == -2)
                    {
                        failed_segment_exists = true;
                        break;
                    }
                    // Also consider the initial_net_repr itself if it was processed as a segment and failed
                    if (initial_net_repr->track_assignment == -2 && gen_net_ptr->parent_representation == initial_net_repr)
                    {
                        failed_segment_exists = true; // Should be covered by initial_net_repr->track_assignment == -2 directly
                        break;
                    }
                }
                if (initial_net_repr->track_assignment == -2)
                    failed_segment_exists = true;

                if (!failed_segment_exists)
                { // If no segment explicitly failed, but it's not marked routed, it means incomplete.
                    std::cerr << "Error/Warning: Net " << initial_net_repr->net_id_str << " was not fully routed. Span ["
                              << initial_net_repr->left_most_x << "," << initial_net_repr->right_most_x
                              << "] is_routed: " << initial_net_repr->is_routed << std::endl;
                }
                all_original_nets_routed = false;
            }
        }

        if (num_tracks_used_ > 0)
            tracks_layout_.resize(num_tracks_used_);
        else
            tracks_layout_.clear();

        std::cout << "Left-Edge Algorithm finished. Tracks used: " << num_tracks_used_ << std::endl;
        return all_original_nets_routed;
    }

    // --- Net Merging Methods ---
    std::vector<NetChannelRepresentation *> DetailedRouter::perform_net_merging(const std::vector<NetChannelRepresentation *> &initial_nets)
    {
        std::cout << "Performing net merging..." << std::endl;
        if (initial_nets.empty())
        {
            return {};
        }
        std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> compatibility_adj;
        build_compatibility_graph(initial_nets, compatibility_adj);
        if (compatibility_adj.empty())
        {
            std::cout << "  No compatible nets found for merging. Proceeding with original nets." << std::endl;
            return initial_nets;
        }
        std::vector<std::vector<const NetChannelRepresentation *>> cliques = find_maximal_cliques(initial_nets, compatibility_adj);
        if (cliques.empty())
        {
            std::cout << "  No cliques found in compatibility graph. Proceeding with original nets." << std::endl;
            return initial_nets;
        }
        std::cout << "  Found " << cliques.size() << " maximal cliques." << std::endl;
        std::vector<NetChannelRepresentation *> final_routable_nets = create_merged_nets_from_cliques(initial_nets, cliques);
        std::cout << "Net merging resulted in " << final_routable_nets.size() << " routable net entities." << std::endl;
        return final_routable_nets;
    }

    void DetailedRouter::build_compatibility_graph(
        const std::vector<NetChannelRepresentation *> &nets,
        std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> &adj) const
    {
        std::cout << "  Building compatibility graph..." << std::endl;
        adj.clear();
        for (size_t i = 0; i < nets.size(); ++i)
        {
            for (size_t j = i + 1; j < nets.size(); ++j)
            {
                const auto *net1 = nets[i];
                const auto *net2 = nets[j];
                if (!net1 || !net2)
                    continue;
                if (net1->left_most_x > net1->right_most_x || net2->left_most_x > net2->right_most_x)
                    continue; // Skip invalid span nets

                bool overlap = std::max(net1->left_most_x, net2->left_most_x) < std::min(net1->right_most_x, net2->right_most_x);
                if (!overlap)
                {
                    adj[net1].insert(net2);
                    adj[net2].insert(net1);
                }
            }
        }
    }

    std::vector<std::vector<const NetChannelRepresentation *>> DetailedRouter::find_maximal_cliques(
        const std::vector<NetChannelRepresentation *> &all_nets,
        const std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> &adj)
    {
        std::cout << "  Finding maximal cliques using Bron-Kerbosch..." << std::endl;
        std::vector<std::vector<const NetChannelRepresentation *>> cliques;
        std::vector<const NetChannelRepresentation *> R;
        std::vector<const NetChannelRepresentation *> P;
        std::vector<const NetChannelRepresentation *> X;

        // Populate P with nodes that are actually in the graph (have entries in adj)
        // This helps focus the algorithm on the relevant part of the graph.
        for (const auto *net : all_nets)
        {
            if (adj.count(net))
            { // Only consider nodes that have compatibility edges
                P.push_back(net);
            }
        }
        // If no nodes have compatibility edges, no cliques can be formed from adj.
        // However, individual nets themselves could be considered cliques of size 1.
        // The current Bron-Kerbosch below finds cliques in the graph defined by `adj`.

        bron_kerbosch_recursive(R, P, X, adj, cliques);
        return cliques;
    }

    void DetailedRouter::bron_kerbosch_recursive(
        std::vector<const NetChannelRepresentation *> &R,
        std::vector<const NetChannelRepresentation *> &P,
        std::vector<const NetChannelRepresentation *> &X,
        const std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> &adj,
        std::vector<std::vector<const NetChannelRepresentation *>> &cliques)
    {

        if (P.empty() && X.empty())
        {
            if (!R.empty())
            {
                cliques.push_back(R);
            }
            return;
        }
        if (P.empty())
            return;

        std::vector<const NetChannelRepresentation *> P_copy = P;
        for (const NetChannelRepresentation *vertex_v : P_copy)
        {
            if (std::find(P.begin(), P.end(), vertex_v) == P.end())
                continue;

            std::vector<const NetChannelRepresentation *> R_new = R;
            R_new.push_back(vertex_v);

            std::vector<const NetChannelRepresentation *> P_new;
            std::vector<const NetChannelRepresentation *> X_new;

            if (adj.count(vertex_v))
            {
                const auto &neighbors_of_v = adj.at(vertex_v);
                for (const NetChannelRepresentation *p_node : P)
                {
                    if (neighbors_of_v.count(p_node))
                    {
                        P_new.push_back(p_node);
                    }
                }
                for (const NetChannelRepresentation *x_node : X)
                {
                    if (neighbors_of_v.count(x_node))
                    {
                        X_new.push_back(x_node);
                    }
                }
            }

            bron_kerbosch_recursive(R_new, P_new, X_new, adj, cliques);

            P.erase(std::remove(P.begin(), P.end(), vertex_v), P.end());
            X.push_back(vertex_v);
        }
    }

    std::vector<NetChannelRepresentation *> DetailedRouter::create_merged_nets_from_cliques(
        const std::vector<NetChannelRepresentation *> &original_nets,
        const std::vector<std::vector<const NetChannelRepresentation *>> &cliques)
    {
        std::cout << "  Creating merged nets from cliques..." << std::endl;
        std::vector<NetChannelRepresentation *> routable_nets_result;
        std::set<const NetChannelRepresentation *> covered_original_nets;

        owned_generated_nets_.clear();
        int merged_net_counter = 0;
        std::vector<std::vector<const NetChannelRepresentation *>> sorted_cliques = cliques;
        std::sort(sorted_cliques.begin(), sorted_cliques.end(), [](const auto &a, const auto &b)
                  { return a.size() > b.size(); });
        for (const auto &clique : sorted_cliques)
        {
            if (clique.size() > 1)
            {
                bool clique_is_useful = false;
                for (const auto *member_net : clique)
                {
                    if (covered_original_nets.find(member_net) == covered_original_nets.end())
                    {
                        clique_is_useful = true;
                        break;
                    }
                }
                if (!clique_is_useful)
                    continue;
                std::string merged_id = "merged_" + std::to_string(merged_net_counter++);
                owned_generated_nets_.emplace_back(std::make_unique<NetChannelRepresentation>(clique, merged_id));
                NetChannelRepresentation *new_merged_net = owned_generated_nets_.back().get();
                routable_nets_result.push_back(new_merged_net);
                for (const auto *member_net : clique)
                {
                    covered_original_nets.insert(member_net);
                }
            }
        }
        for (NetChannelRepresentation *orig_net : original_nets)
        {
            if (covered_original_nets.find(orig_net) == covered_original_nets.end())
            {
                routable_nets_result.push_back(orig_net);
            }
        }
        if (routable_nets_result.empty() && !original_nets.empty())
        {
            return original_nets;
        }
        return routable_nets_result;
    }

    // --- Dogleg Helper Stub ---
    NetChannelRepresentation *DetailedRouter::attempt_dogleg_split(
        NetChannelRepresentation *current_segment_orig,
        int proposed_dogleg_x,
        int track_idx_for_first_part)
    {
        std::cout << "Attempting dogleg split for net: " << current_segment_orig->net_id_str
                  << " at x=" << proposed_dogleg_x << " on track " << track_idx_for_first_part << std::endl;
        // 1. Validate proposed_dogleg_x: must be > current_segment_orig->current_segment_start_x and < current_segment_orig->current_segment_end_x
        //    and ideally > any pin x_coord in the first part.
        if (!current_segment_orig ||
            proposed_dogleg_x <= current_segment_orig->current_segment_start_x ||
            proposed_dogleg_x >= current_segment_orig->current_segment_end_x)
        {
            // std::cerr << "  Invalid proposed_dogleg_x: " << proposed_dogleg_x << " for segment span ["
            //           << current_segment_orig->current_segment_start_x << ", " << current_segment_orig->current_segment_end_x << "]" << std::endl;
            return nullptr;
        }

        // 2. Try to place the first part of the segment [current_start, proposed_dogleg_x] on track_idx_for_first_part
        //    We need to use assign_net_to_track or similar logic. assign_net_to_track modifies the net_repr, so be careful.
        //    For now, assume if we are here, the track is viable for *some* part.
        //    A simpler check for this stub: can_place for the first segment.
        if (!tracks_layout_[track_idx_for_first_part].can_place(current_segment_orig, current_segment_orig->current_segment_start_x, proposed_dogleg_x))
        {
            // std::cout << "  Cannot place first part of dogleg for " << current_segment_orig->net_id_str << " on track " << track_idx_for_first_part << std::endl;
            return nullptr;
        }
        tracks_layout_[track_idx_for_first_part].add_segment(current_segment_orig, current_segment_orig->current_segment_start_x, proposed_dogleg_x);
        current_segment_orig->track_assignment = track_idx_for_first_part; // First part is on this track
        // Mark this segment as partially routed up to proposed_dogleg_x
        // The is_routed flag on current_segment_orig might mean fully routed. We need a more nuanced state or rely on segment end.
        // Let's update its current_segment_end_x. The original right_most_x remains unchanged.
        int original_first_part_end_x = current_segment_orig->current_segment_end_x; // Save old end
        current_segment_orig->current_segment_end_x = proposed_dogleg_x;
        // current_segment_orig->is_routed = false; // It's not fully routed yet if it needs a dogleg

        // 3. Create new NetChannelRepresentation for the second part [proposed_dogleg_x, original_full_end_x]
        //    The new segment effectively starts *after* the jog at proposed_dogleg_x.
        //    For routing, its pins are those >= proposed_dogleg_x from the original net.
        std::vector<ChannelPin> second_segment_pins;
        NetChannelRepresentation *ultimate_parent = current_segment_orig->parent_representation ? current_segment_orig->parent_representation : current_segment_orig;
        for (const auto &pin : ultimate_parent->pins)
        {
            if (pin.x_coord >= proposed_dogleg_x)
            { // Pins for the second segment
                second_segment_pins.push_back(pin);
            }
        }

        if (second_segment_pins.empty() && proposed_dogleg_x < ultimate_parent->right_most_x)
        {
            // This can happen if the split point is after the last pin but before the net's actual end.
            // The second segment still needs to span to right_most_x.
            // The new segment might have no pins of its own but needs to exist for connectivity.
            // For simplicity, a segment should have at least one connection point or cover a span with pins.
            // If no pins remain for the second segment, it means the net essentially ended at or before proposed_dogleg_x from a pin perspective.
            // Or it means the remaining part is just a horizontal track to the rightmost extent without further pins.
            // Consider if such a segment is valid or if the original segment should just end at proposed_dogleg_x.
            // For now, if no pins for second segment but original net extends further, this is tricky.
            // Let's assume for a valid dogleg, the second part should connect to something or complete the span.
            // If second_segment_pins is empty but proposed_dogleg_x < original_first_part_end_x, it is an issue.
            std::cout << "  Warning: Second part of dogleg for " << current_segment_orig->net_id_str << " has no pins but original net extends to " << ultimate_parent->right_most_x << std::endl;
            // Potentially, this segment doesn't need explicit pins if it's just a continuation.
            // Its start is proposed_dogleg_x, end is ultimate_parent->right_most_x.
        }

        // The new segment starts conceptually at proposed_dogleg_x (or column after jog)
        // and ends at the original net's true right_most_x.
        std::string new_segment_id = ultimate_parent->net_id_str + "_dogleg_" + std::to_string(owned_generated_nets_.size());
        owned_generated_nets_.emplace_back(std::make_unique<NetChannelRepresentation>(
            ultimate_parent, // The new segment is a child of the original (or ultimate parent of current segment)
            new_segment_id,
            proposed_dogleg_x,             // Second segment starts at dogleg point
            ultimate_parent->right_most_x, // And goes to the original full end
            second_segment_pins));
        NetChannelRepresentation *second_segment = owned_generated_nets_.back().get();
        second_segment->is_primary_segment = false; // Clearly a subsequent segment

        std::cout << "  Created dogleg segment: " << second_segment->net_id_str << " span ["
                  << second_segment->current_segment_start_x << ", " << second_segment->current_segment_end_x << "]" << std::endl;

        // The original segment (current_segment_orig) is now shorter.
        // If current_segment_orig->current_segment_end_x == current_segment_orig->right_most_x (based on its pins),
        // then it might be considered fully routed for its part.
        // The overall net (ultimate_parent) is routed when all its pins are covered by routed segments.
        if (current_segment_orig->current_segment_end_x >= current_segment_orig->right_most_x)
        {
            current_segment_orig->is_routed = true; // This segment is now fully routed.
        }

        return second_segment; // Return the new segment to be added to routing queue
    }

    bool DetailedRouter::route_channel(Channel &channel_data)
    {
        if (channel_data.nets.empty())
        {
            std::cout << "Channel " << channel_data.id << ": No nets to route." << std::endl;
            return true;
        }

        channel_length_ = channel_data.width;
        channel_capacity_ = channel_data.height;
        num_tracks_used_ = 0;
        tracks_layout_.clear();
        vcg_above_.clear();
        vcg_below_.clear();
        owned_generated_nets_.clear();

        std::cout << "Starting detailed routing for channel: " << channel_data.id
                  << " Length: " << channel_length_ << " Capacity: " << channel_capacity_
                  << " Original Nets: " << channel_data.nets.size() << std::endl;

        initialize_net_representations(channel_data.nets);

        std::sort(channel_data.nets.begin(), channel_data.nets.end(), [](const NetChannelRepresentation *a, const NetChannelRepresentation *b)
                  {
            if (a->left_most_x != b->left_most_x) {
                return a->left_most_x < b->left_most_x;
            }
            return a->net_id_str < b->net_id_str; });

        std::vector<NetChannelRepresentation *> routable_nets = perform_net_merging(channel_data.nets);

        // If merging produces an empty list but original nets existed, fall back to original nets.
        // This ensures LEA always has something to work on if there were initial nets.
        if (routable_nets.empty() && !channel_data.nets.empty())
        {
            std::cout << "Warning: Net merging resulted in an empty set of routable nets. Falling back to original nets." << std::endl;
            routable_nets = channel_data.nets; // Use original nets if merging produced nothing (e.g. all nets filtered)
        }

        std::sort(routable_nets.begin(), routable_nets.end(), [](const NetChannelRepresentation *a, const NetChannelRepresentation *b)
                  {
            if (a->left_most_x != b->left_most_x) {
                return a->left_most_x < b->left_most_x;
            }
            if (a->right_most_x != b->right_most_x) {
                return a->right_most_x < b->right_most_x;
            }
            return a->net_id_str < b->net_id_str; });

        std::cout << "Proceeding to VCG and LEA with " << routable_nets.size() << " routable entities." << std::endl;
        if (routable_nets.empty())
        { // Double check after potential fallback
            std::cout << "No routable net entities to process after merging. Skipping VCG and LEA." << std::endl;
            // Depending on desired behavior, could return true if original nets were also empty,
            // or false if original nets existed but couldn't be processed.
            return channel_data.nets.empty();
        }

        build_vertical_constraint_graph(routable_nets);

        bool lea_success = left_edge_algorithm_route(routable_nets, channel_data);
        if (!lea_success)
        {
            std::cerr << "Channel routing failed for channel " << channel_data.id << " using LEA." << std::endl;
            return false;
        }

        std::cout << "Channel " << channel_data.id << " routed. Tracks used: " << num_tracks_used_ << std::endl;
        return true;
    }

} // namespace router