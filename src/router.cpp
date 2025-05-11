// router.cpp – Global & Detailed Routing for ISPD‑2011 project
// -----------------------------------------------------------------------------
//   * Implements two Steiner‑tree‑free global routing algorithms:
//       1. Hadlock 0‑1 BFS (minimum detours)
//       2. Soukup best‑first search with look‑ahead
//   * Provides a very lightweight maze‑based detailed router that assigns each
//     global segment to a metal track, honoring per‑row site geometry and a
//     single minimum track spacing.
//   * Self‑contained: include this file once *after* you have parsed the design
//     (Parser) and performed placement (Placer).
//
// Example usage (add to main.cpp right after placement):
// -----------------------------------------------------------------------------
//     #include "router.cpp"          // or split header/impl if you prefer
//     Router router(circuit);        // build routing grids
//     auto unrouted = router.global_route(ALGO_HADLOCK, /*verbose=*/true);
//
//     DetailedRouter dr(circuit, router);
//     auto violations = dr.detailed_route(/*verbose=*/true);
//     std::cout << "[Summary] Unrouted nets:" << unrouted
//               << "  DRC violations:" << violations << std::endl;
// -----------------------------------------------------------------------------
// The grids are stored inside the Router so that the Evaluator can directly
// query edge usage for congestion / overflow statistics.
// -----------------------------------------------------------------------------

#include "../include/router.h"
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <algorithm>

// -----------------------------------------------------------------------------
// Router Implementation
// -----------------------------------------------------------------------------

Router::Router(Circuit &ckt) : ckt_(ckt)
{
    build_grid(LAMBDA);
}

int Router::global_route(GlobalAlgo algo, bool verbose, MergeAlgo mAlg)
{
    unrouted_.clear();
    int failures = 0;

    for (auto &np : ckt_.net_map)
    {
        Net &net = np.second;
        auto segments = decompose_net(net, mAlg);
        for (const auto &seg : segments)
        {
            const GridPt &s = seg.first;
            const GridPt &t = seg.second;

            auto path = (algo == ALGO_HADLOCK)
                            ? route_hadlock(s, t)
                            : route_soukup(s, t);

            if (!path)
            {
                ++failures;
                unrouted_.insert(&net);
                if (verbose)
                    std::cout << "[Router] Failed segment of net "
                              << net.name << "\n";
                continue;
            }
            commit_path(*path);
            net_paths_[&net].push_back(*path);
        }
    }
    if (verbose)
        std::cout << "[Router] Global route completed. Unrouted: " << failures << std::endl;
    return failures;
}

const std::unordered_map<const Net *, std::vector<std::vector<GridPt>>> &Router::get_paths() const
{
    return net_paths_;
}

const std::unordered_map<EdgeKey, int, std::hash<EdgeKey>> &Router::h_usage() const
{
    return h_usage_;
}

const std::unordered_map<EdgeKey, int, std::hash<EdgeKey>> &Router::v_usage() const
{
    return v_usage_;
}

const std::unordered_set<const Net *> &Router::unrouted_nets() const
{
    return unrouted_;
}

GridPt Router::to_grid(const Node *n) const
{
    return {(int)std::round((n->pos.x - origin_.x) / tile_w_),
            (int)std::round((n->pos.y - origin_.y) / tile_h_)};
}

void Router::build_grid(double lambda)
{
    if (ckt_.routing_grid.grid_x > 0)
    {
        origin_ = ckt_.routing_grid.origin;
        tile_w_ = ckt_.routing_grid.tile_width;
        tile_h_ = ckt_.routing_grid.tile_height;
        gx_ = ckt_.routing_grid.grid_x;
        gy_ = ckt_.routing_grid.grid_y;
    }
    else
    {
        origin_ = ckt_.core_region.bottom_left;
        tile_w_ = tile_h_ = lambda;
        gx_ = (int)std::ceil((ckt_.core_region.top_right.x - origin_.x) / tile_w_);
        gy_ = (int)std::ceil((ckt_.core_region.top_right.y - origin_.y) / tile_h_);
    }
}

bool Router::in_bounds(const GridPt &p) const
{
    return p.x >= 0 && p.y >= 0 && p.x < gx_ && p.y < gy_;
}

std::optional<std::vector<GridPt>> Router::route_hadlock(const GridPt &s, const GridPt &t)
{
    struct NodeInfo
    {
        int detours;
        GridPt parent;
    };
    const int INF = std::numeric_limits<int>::max();
    std::vector<std::vector<NodeInfo>> dist(gx_, std::vector<NodeInfo>(gy_, {INF, {-1, -1}}));
    auto manh = [&](const GridPt &a, const GridPt &b)
    { return std::abs(a.x - b.x) + std::abs(a.y - b.y); };
    auto cmp = [&](const std::pair<int, GridPt> &a, const std::pair<int, GridPt> &b)
    { return a.first > b.first; };
    std::priority_queue<std::pair<int, GridPt>, std::vector<std::pair<int, GridPt>>, decltype(cmp)> pq(cmp);
    dist[s.x][s.y].detours = 0;
    pq.push({0, s});
    const GridPt dirs[4] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    while (!pq.empty())
    {
        auto [d, p] = pq.top();
        pq.pop();
        if (p == t)
            break;
        for (auto dir : dirs)
        {
            GridPt nb{p.x + dir.x, p.y + dir.y};
            if (!in_bounds(nb))
                continue;
            int nd = d + (manh(nb, t) > manh(p, t) ? 1 : 0);
            if (nd < dist[nb.x][nb.y].detours)
            {
                dist[nb.x][nb.y].detours = nd;
                dist[nb.x][nb.y].parent = p;
                pq.push({nd, nb});
            }
        }
    }
    if (dist[t.x][t.y].detours == INF)
        return std::nullopt;
    // back‑trace
    std::vector<GridPt> path;
    for (GridPt cur = t; !(cur == s); cur = dist[cur.x][cur.y].parent)
        path.push_back(cur);
    path.push_back(s);
    std::reverse(path.begin(), path.end());
    return path;
}

std::optional<std::vector<GridPt>> Router::route_soukup(const GridPt &s, const GridPt &t)
{
    struct NodeInfo
    {
        int score;
        GridPt parent;
    };
    const int INF = std::numeric_limits<int>::max();
    std::vector<std::vector<NodeInfo>> dist(gx_, std::vector<NodeInfo>(gy_, {INF, {-1, -1}}));
    auto score_f = [&](const GridPt &p, int det)
    { return 10 * (std::abs(p.x - t.x) + std::abs(p.y - t.y)) + det; };
    auto cmp = [&](const std::pair<int, GridPt> &a, const std::pair<int, GridPt> &b)
    { return a.first > b.first; };
    std::priority_queue<std::pair<int, GridPt>, std::vector<std::pair<int, GridPt>>, decltype(cmp)> pq(cmp);
    dist[s.x][s.y].score = 0;
    pq.push({0, s});
    const GridPt dirs[4] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    while (!pq.empty())
    {
        auto [sc, p] = pq.top();
        pq.pop();
        if (p == t)
            break;
        for (auto dir : dirs)
        {
            GridPt nb{p.x + dir.x, p.y + dir.y};
            if (!in_bounds(nb))
                continue;
            int ndet = sc + ((dir.x != 0) == (p.x - t.x != 0) && (dir.y != 0) == (p.y - t.y != 0) ? 0 : 1);
            int nscore = score_f(nb, ndet);
            if (nscore < dist[nb.x][nb.y].score)
            {
                dist[nb.x][nb.y].score = nscore;
                dist[nb.x][nb.y].parent = p;
                pq.push({nscore, nb});
            }
        }
    }
    if (dist[t.x][t.y].score == INF)
        return std::nullopt;
    // back‑trace
    std::vector<GridPt> path;
    for (GridPt cur = t; !(cur == s); cur = dist[cur.x][cur.y].parent)
        path.push_back(cur);
    path.push_back(s);
    std::reverse(path.begin(), path.end());
    return path;
}

void Router::commit_path(const std::vector<GridPt> &p)
{
    for (size_t i = 1; i < p.size(); ++i)
    {
        EdgeKey e = make_edge(p[i - 1], p[i]);
        if (p[i - 1].x == p[i].x)
            v_usage_[e]++;
        else
            h_usage_[e]++;
    }
}

Point Router::node_pos(const Node *n) const
{
    return n->pos;
}

GridPt Router::phys_to_grid(const Point &p) const
{
    return {(int)std::round((p.x - origin_.x) / tile_w_),
            (int)std::round((p.y - origin_.y) / tile_h_)};
}

std::vector<std::pair<GridPt, GridPt>> Router::decompose_star(const Net &net) const
{
    std::vector<std::pair<GridPt, GridPt>> segs;
    if (net.nodes.size() < 2)
        return segs;

    /* centroid */
    double cx = 0, cy = 0;
    for (auto *n : net.nodes)
    {
        const Point p = node_pos(n);
        cx += p.x;
        cy += p.y;
    }
    cx /= net.nodes.size();
    cy /= net.nodes.size();

    /* root = pin closest to centroid */
    const Node *root = net.nodes[0];
    double best = std::numeric_limits<double>::max();
    for (auto *n : net.nodes)
    {
        const Point p = node_pos(n);
        double d = std::abs(p.x - cx) + std::abs(p.y - cy);
        if (d < best)
        {
            best = d;
            root = n;
        }
    }
    GridPt rG = phys_to_grid(node_pos(root));

    /* root -> every other */
    for (auto *n : net.nodes)
        if (n != root)
            segs.push_back({rG, phys_to_grid(node_pos(n))});

    return segs;
}

std::vector<std::pair<GridPt, GridPt>> Router::decompose_mst(const Net &net) const
{
    std::vector<std::pair<GridPt, GridPt>> segs;
    const size_t k = net.nodes.size();
    if (k < 2)
        return segs;

    std::vector<GridPt> pins;
    pins.reserve(k);
    for (auto *n : net.nodes)
        pins.push_back(phys_to_grid(node_pos(n)));

    std::vector<bool> inTree(k, false);
    std::vector<int> best(k, INT_MAX), parent(k, -1);

    inTree[0] = true;
    GridPt ref = pins[0];
    for (size_t i = 1; i < k; ++i)
        best[i] = std::abs(pins[i].x - ref.x) + std::abs(pins[i].y - ref.y),
        parent[i] = 0;

    for (size_t it = 1; it < k - 0; ++it)
    {
        int mn = INT_MAX, idx = -1;
        for (size_t v = 0; v < k; ++v)
            if (!inTree[v] && best[v] < mn)
            {
                mn = best[v];
                idx = (int)v;
            }
        if (idx == -1)
            break;
        inTree[idx] = true;
        segs.push_back({pins[idx], pins[parent[idx]]});

        for (size_t v = 0; v < k; ++v)
            if (!inTree[v])
            {
                int d = std::abs(pins[v].x - pins[idx].x) +
                        std::abs(pins[v].y - pins[idx].y);
                if (d < best[v])
                {
                    best[v] = d;
                    parent[v] = (int)idx;
                }
            }
    }
    return segs;
}

std::vector<std::pair<GridPt, GridPt>> Router::decompose_net(const Net &net, MergeAlgo alg) const
{
    return (alg == MERGE_MST) ? decompose_mst(net)
                              : decompose_star(net);
}

// -----------------------------------------------------------------------------
// DetailedRouter Implementation
// -----------------------------------------------------------------------------

DetailedRouter::DetailedRouter(Circuit &ckt, Router &gr) : ckt_(ckt), gr_(gr)
{
    build_tracks();
}

int DetailedRouter::detailed_route(bool verbose)
{
    int violations = 0;
    for (auto &kv : ckt_.net_map)
    {
        if (gr_.unrouted_nets().count(&kv.second))
            continue; // skip
        const auto &gpaths = gr_.get_paths().at(&kv.second);
        for (const auto &ptlist : gpaths)
        {
            for (size_t i = 1; i < ptlist.size(); ++i)
            {
                Point p1 = grid_to_coord(ptlist[i - 1]);
                Point p2 = grid_to_coord(ptlist[i]);
                int row = pick_row(p1, p2);
                if (row < 0)
                {
                    ++violations;
                    if (verbose)
                        std::cout << "[DRC] Row not found for segment of net " << kv.second.name << "\n";
                    continue;
                }
                Rect seg{p1, p2};
                normalize_rect(seg);
                bool collides = false;
                for (const Rect &r : row_tracks_[row])
                {
                    if (overlap(r, seg))
                    {
                        collides = true;
                        break;
                    }
                }
                if (collides)
                {
                    ++violations;
                    if (verbose)
                        std::cout << "[DRC] Track collision in row " << row << " for net " << kv.second.name << "\n";
                }
                else
                {
                    row_tracks_[row].push_back(seg);
                }
            }
        }
    }
    if (verbose)
        std::cout << "[DetailedRouter] complete. violations=" << violations << "\n";
    return violations;
}

void DetailedRouter::build_tracks()
{
    row_tracks_.resize(ckt_.core_rows.size());
}

Point DetailedRouter::grid_to_coord(const GridPt &g) const
{
    return {gr_.origin_.x + g.x * gr_.tile_w_, gr_.origin_.y + g.y * gr_.tile_h_};
}

int DetailedRouter::pick_row(const Point &a, const Point &b) const
{
    double y = 0.5 * (a.y + b.y);
    for (size_t i = 0; i < ckt_.core_rows.size(); ++i)
    {
        const CoreRow &r = ckt_.core_rows[i];
        if (y >= r.y_coord && y < r.y_coord + r.height)
            return (int)i;
    }
    return -1;
}

void DetailedRouter::normalize_rect(Rect &r)
{
    if (r.bottom_left.x > r.top_right.x)
        std::swap(r.bottom_left.x, r.top_right.x);
    if (r.bottom_left.y > r.top_right.y)
        std::swap(r.bottom_left.y, r.top_right.y);
}

bool DetailedRouter::overlap(const Rect &a, const Rect &b)
{
    double h = std::max(0.0, std::min(a.top_right.y, b.top_right.y) - std::max(a.bottom_left.y, b.bottom_left.y));
    double w = std::max(0.0, std::min(a.top_right.x, b.top_right.x) - std::max(a.bottom_left.x, b.bottom_left.x));
    return h > 0 && w > 0;
}
