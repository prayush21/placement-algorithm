#ifndef ROUTER_H
#define ROUTER_H

#include "structures.h"
#include <queue>
#include <unordered_map>
#include <vector>
#include <deque>
#include <limits>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <optional>
#include <unordered_set>

// Helper data structure definitions
struct GridPt
{
    int x, y;
    bool operator==(const GridPt &o) const { return x == o.x && y == o.y; }
};

namespace std
{
    template <>
    struct hash<GridPt>
    {
        size_t operator()(const GridPt &p) const noexcept
        {
            return hash<int>()(p.x) ^ (hash<int>()(p.y) << 1);
        }
    };

    template <>
    struct hash<pair<GridPt, GridPt>>
    {
        size_t operator()(const pair<GridPt, GridPt> &p) const noexcept
        {
            size_t h1 = hash<GridPt>()(p.first);
            size_t h2 = hash<GridPt>()(p.second);
            return h1 ^ (h2 << 1);
        }
    };
}

// Forward declarations and type aliases
using EdgeKey = std::pair<GridPt, GridPt>;

static inline EdgeKey make_edge(const GridPt &a, const GridPt &b)
{
    return (a.x < b.x || (a.x == b.x && a.y < b.y)) ? EdgeKey{a, b} : EdgeKey{b, a};
}

enum GlobalAlgo
{
    ALGO_HADLOCK,
    ALGO_SOUKUP
};

enum MergeAlgo
{
    MERGE_STAR,
    MERGE_MST
};

class Router
{
public:
    explicit Router(Circuit &ckt);
    int global_route(GlobalAlgo algo = ALGO_HADLOCK, bool verbose = false, MergeAlgo mAlg = MERGE_STAR);
    const std::unordered_map<const Net *, std::vector<std::vector<GridPt>>> &get_paths() const;
    const std::unordered_map<EdgeKey, int, std::hash<EdgeKey>> &h_usage() const;
    const std::unordered_map<EdgeKey, int, std::hash<EdgeKey>> &v_usage() const;
    const std::unordered_set<const Net *> &unrouted_nets() const;

    void save_global_routing_results(const std::string &filename) const;
    std::optional<int> load_global_routing_results(const std::string &filename);

    friend class DetailedRouter;

protected:
    Circuit &ckt_;
    Point origin_;
    double tile_w_ = 1.0, tile_h_ = 1.0;
    int gx_ = 0, gy_ = 0;

    std::unordered_map<const Net *, std::vector<std::vector<GridPt>>> net_paths_;
    std::unordered_map<EdgeKey, int, std::hash<EdgeKey>> h_usage_;
    std::unordered_map<EdgeKey, int, std::hash<EdgeKey>> v_usage_;
    std::unordered_set<const Net *> unrouted_;

    void build_grid(double lambda);
    bool in_bounds(const GridPt &p) const;
    std::optional<std::vector<GridPt>> route_hadlock(const GridPt &s, const GridPt &t);
    std::optional<std::vector<GridPt>> route_soukup(const GridPt &s, const GridPt &t);
    void commit_path(const std::vector<GridPt> &p);
    Point node_pos(const Node *n) const;
    GridPt to_grid(const Node *n) const;
    GridPt phys_to_grid(const Point &p) const;
    std::vector<std::pair<GridPt, GridPt>> decompose_star(const Net &net) const;
    std::vector<std::pair<GridPt, GridPt>> decompose_mst(const Net &net) const;
    std::vector<std::pair<GridPt, GridPt>> decompose_net(const Net &net, MergeAlgo alg) const;
};

class DetailedRouter
{
public:
    DetailedRouter(Circuit &ckt, Router &gr);
    int detailed_route(bool verbose = false);

private:
    Circuit &ckt_;
    Router &gr_;
    std::vector<std::vector<Rect>> row_tracks_;
    const double MIN_SPACING = 1.0 * LAMBDA;

    void build_tracks();
    Point grid_to_coord(const GridPt &g) const;
    int pick_row(const Point &a, const Point &b) const;
    void normalize_rect(Rect &r);
    bool overlap(const Rect &a, const Rect &b);
};

#endif // ROUTER_H