#ifndef DETAILED_ROUTER_H
#define DETAILED_ROUTER_H

#include "structures.h" // For Channel, Net, PinSide, ChannelPin, NetChannelRepresentation
#include <vector>
#include <string>
#include <map>       // For VCG
#include <set>       // For VCG
#include <algorithm> // For std::sort, std::min, std::max
#include <memory>    // For std::unique_ptr
// #include <iostream> // For std::cout, std::cerr - remove if only for .cpp

namespace router
{

    // PinSide, ChannelPin, and NetChannelRepresentation are now defined in structures.h
    // and included via "structures.h".
    // Remove their definitions from here.

    class DetailedRouter
    {
    public:
        DetailedRouter();
        ~DetailedRouter();

        // Main routing function expects Channel to be populated with NetChannelRepresentation objects
        bool route_channel(Channel &channel_data);

    private:
        // Internal representation of a net segment being placed on a track
        struct TrackSegment
        {
            NetChannelRepresentation *net_repr; // Pointer to the net representation from Channel.nets
            int x_start;                        // Start x-coordinate of this segment on the track
            int x_end;                          // End x-coordinate of this segment on the track

            TrackSegment(NetChannelRepresentation *repr, int start, int end) : net_repr(repr), x_start(start), x_end(end) {}
        };

        // Represents a single track in the channel
        struct Track
        {
            std::vector<TrackSegment> segments; // Segments placed on this track
            int track_idx;                      // Its own index
            Track() : track_idx(0) {}           // Default constructor

            Track(int idx) : track_idx(idx) {}

            // Check if a new segment can be added to this track
            bool can_place(const NetChannelRepresentation *net_to_place, int seg_start_x, int seg_end_x) const
            {
                // Check for overlaps with existing segments on this track
                for (const auto &existing_segment : segments)
                {
                    // Check for horizontal overlap: max(start1, start2) < min(end1, end2)
                    if (std::max(existing_segment.x_start, seg_start_x) < std::min(existing_segment.x_end, seg_end_x))
                    {
                        return false; // Overlap detected
                    }
                }
                return true; // No overlap, can place
            }

            void add_segment(NetChannelRepresentation *net_to_add, int seg_start_x, int seg_end_x)
            {
                segments.emplace_back(net_to_add, seg_start_x, seg_end_x);
                // Optional: Keep segments sorted by x_start for easier analysis or merging later
                std::sort(segments.begin(), segments.end(), [](const TrackSegment &a, const TrackSegment &b)
                          { return a.x_start < b.x_start; });
            }
        };

        std::vector<Track> tracks_layout_; // Stores the assignment of nets to tracks
        int channel_length_;               // Length of the channel (width from .chn)
        int channel_capacity_;             // Max tracks available (height from .chn, or dynamically determined)
        int num_tracks_used_;

        // Processed nets for routing. This might be a copy or direct pointers from channel_data.nets
        // If channel_data.nets are modified directly (e.g. track_assignment), this might not be needed
        // or could store supplemental router-specific state.
        // For now, assume we operate on channel_data.nets[i] directly for routing state like track_assignment.
        // std::vector<NetChannelRepresentation> processed_nets_; // Let's try to work with channel_data.nets directly

        // VCG: An edge U -> V in VCG means U must be in a track above V.
        // Key: NetRepresentation that is BELOW, Value: Set of NetRepresentations that must be ABOVE.
        std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> vcg_above_; // Stores nets that must be above the key net
        std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> vcg_below_; // Stores nets that must be below the key net
        std::map<const NetChannelRepresentation *, int> vcg_in_degree_;                                    // For LEA processing: number of constraints on this net FROM UNROUTED NETS

        // Renamed member to own NetChannelRepresentation objects created for merged nets AND dogleg segments
        std::vector<std::unique_ptr<NetChannelRepresentation>> owned_generated_nets_;

        // 1. Initialization
        void initialize_net_representations(const std::vector<NetChannelRepresentation *> &nets); // Operates on a list of nets

        // 2. Net Merging (Preprocessing)
        std::vector<NetChannelRepresentation *> perform_net_merging(const std::vector<NetChannelRepresentation *> &initial_nets);
        void build_compatibility_graph(const std::vector<NetChannelRepresentation *> &nets,
                                       std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> &adj) const;
        std::vector<std::vector<const NetChannelRepresentation *>> find_maximal_cliques(
            const std::vector<NetChannelRepresentation *> &all_nets, // Need all nets for iteration context in BronKerbosch
            const std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> &adj);
        void bron_kerbosch_recursive(std::vector<const NetChannelRepresentation *> &R,
                                     std::vector<const NetChannelRepresentation *> &P,
                                     std::vector<const NetChannelRepresentation *> &X,
                                     const std::map<const NetChannelRepresentation *, std::set<const NetChannelRepresentation *>> &adj,
                                     std::vector<std::vector<const NetChannelRepresentation *>> &cliques);
        std::vector<NetChannelRepresentation *> create_merged_nets_from_cliques(
            const std::vector<NetChannelRepresentation *> &original_nets,
            const std::vector<std::vector<const NetChannelRepresentation *>> &cliques);

        // 3. Vertical Constraint Graph (VCG)
        void build_vertical_constraint_graph(const std::vector<NetChannelRepresentation *> &routable_nets);
        // void calculate_vcg_in_degrees(); // May not be needed if LEA handles it dynamically

        // 4. Modified Left-Edge Algorithm (LEA) with Doglegs
        bool left_edge_algorithm_route(const std::vector<NetChannelRepresentation *> &routable_nets, Channel &channel_data_ref_for_length_capacity);

        // Helper for LEA: find next net to route
        NetChannelRepresentation *select_next_net_for_lea(std::vector<NetChannelRepresentation *> &active_nets); // pass by ref to modify

        // Helper for LEA: assign net (or segment) to a track
        bool assign_net_to_track(NetChannelRepresentation *net_repr, int track_idx, int segment_start_x, int segment_end_x);

        // Helper for dogleg creation
        // Tries to place the first part of current_segment up to proposed_dogleg_x on track_for_first_part.
        // If successful, updates current_segment, creates a new NetChannelRepresentation for the remainder, adds it to owned_generated_nets_,
        // and returns a pointer to the new (second) segment. Otherwise, returns nullptr.
        NetChannelRepresentation *attempt_dogleg_split(
            NetChannelRepresentation *current_segment_orig, // The segment to be split
            int proposed_dogleg_x,                          // Column *before* which the first part ends / jog occurs
            int track_idx_for_first_part);

        // 5. Net Merging (Placeholder - typically post-LEA or integrated)
        // The user's request implies pre-LEA merging. This can be removed or re-purposed if needed.
        // void optimize_by_merging();
    };

} // namespace router

#endif // DETAILED_ROUTER_H