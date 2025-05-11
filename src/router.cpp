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
#include <fstream>
#include <sstream>

// -----------------------------------------------------------------------------
// Router Implementation
// -----------------------------------------------------------------------------

Router::Router(Circuit &ckt) : ckt_(ckt)
{
    // Log the value of LAMBDA being passed to build_grid
    std::cout << "[Router-Debug] Router constructor: LAMBDA constant value = " << LAMBDA << std::endl;
    build_grid(LAMBDA);
}

int Router::global_route(GlobalAlgo algo, bool verbose, MergeAlgo mAlg)
{
    if (verbose)
        std::cout << "[Router] Starting global_route method." << std::endl;

    unrouted_.clear();
    int failures = 0;

    if (verbose)
        std::cout << "[Router] Processing " << ckt_.net_map.size() << " nets." << std::endl;

    for (auto &np : ckt_.net_map)
    {
        Net &net = np.second;
        if (verbose)
            std::cout << "[Router] Processing net: " << net.name << std::endl;

        if (verbose && net.name == "n468917")
        { // Log specific net's pin coords
            std::cout << "[Router] Physical pin coordinates for net " << net.name << ":" << std::endl;
            for (const auto *node : net.nodes)
            {
                if (node)
                {
                    std::cout << "  Pin: (" << node->pos.x << ", " << node->pos.y << ")" << std::endl;
                }
                else
                {
                    std::cout << "  Pin: (null node)" << std::endl;
                }
            }
        }

        if (verbose)
            std::cout << "[Router] Decomposing net: " << net.name << " using algorithm: " << (mAlg == MERGE_MST ? "MST" : "STAR") << std::endl;
        auto segments = decompose_net(net, mAlg);
        if (verbose)
            std::cout << "[Router] Net " << net.name << " decomposed into " << segments.size() << " segments." << std::endl;

        for (const auto &seg : segments)
        {
            const GridPt &s = seg.first;
            const GridPt &t = seg.second;

            if (verbose)
                std::cout << "[Router] Routing segment for net " << net.name << " from (" << s.x << "," << s.y << ") to (" << t.x << "," << t.y << ")" << std::endl;

            auto path = (algo == ALGO_HADLOCK)
                            ? route_hadlock(s, t)
                            : route_soukup(s, t);

            if (!path)
            {
                ++failures;
                unrouted_.insert(&net);
                if (verbose)
                    std::cout << "[Router] Failed segment of net "
                              << net.name << " from (" << s.x << "," << s.y << ") to (" << t.x << "," << t.y << ")" << "";
                continue;
            }
            if (verbose)
                std::cout << "[Router] Successfully routed segment for net " << net.name << ". Path length: " << path->size() << std::endl;
            commit_path(*path);
            if (verbose)
                std::cout << "[Router] Committed path for segment of net " << net.name << std::endl;
            net_paths_[&net].push_back(*path);
        }
        if (verbose)
            std::cout << "[Router] Finished processing net: " << net.name << std::endl;
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

void Router::build_grid(double lambda_param)
{
    std::cout << "[Router-Debug] build_grid called with lambda_param = " << lambda_param << std::endl;

    if (ckt_.routing_grid.grid_x > 0 && ckt_.routing_grid.grid_y > 0) // Ensure both Gx and Gy from .route are valid
    {
        origin_ = ckt_.routing_grid.origin;
        tile_w_ = ckt_.routing_grid.tile_width;
        tile_h_ = ckt_.routing_grid.tile_height;
        gx_ = ckt_.routing_grid.grid_x;
        gy_ = ckt_.routing_grid.grid_y;

        std::cout << "[Router-Debug] build_grid (initial from .route data):\n"
                  << "  Origin: (" << origin_.x << ", " << origin_.y << ")\n"
                  << "  Tile W/H (from ckt_.routing_grid): (" << tile_w_ << ", " << tile_h_ << ")\n"
                  << "  Grid Dims (gx, gy): (" << gx_ << ", " << gy_ << ")" << std::endl;

        // If tile dimensions from .route are invalid (0 or negative), try to calculate them
        // using core region dimensions and Gx/Gy from .route data.
        if (tile_w_ <= 0 || tile_h_ <= 0)
        {
            std::cout << "[Router-Warning] Tile W/H from .route file are invalid or zero (" << tile_w_ << ", " << tile_h_ << ").\n"
                      << "                 Attempting to calculate from core region and Gx/Gy." << std::endl;

            // Ensure core region dimensions are sensible and Gx/Gy are positive for calculation
            double core_width = ckt_.core_region.top_right.x - origin_.x;  // Width relative to grid origin
            double core_height = ckt_.core_region.top_right.y - origin_.y; // Height relative to grid origin

            if (gx_ > 0 && gy_ > 0 && core_width > 0 && core_height > 0)
            {
                tile_w_ = core_width / gx_;
                tile_h_ = core_height / gy_;
                std::cout << "[Router-Debug] build_grid (recalculated Tile W/H for .route case):\n"
                          << "  Core Width (rel to origin): " << core_width << " (TR.x " << ckt_.core_region.top_right.x << " - Org.x " << origin_.x << ")\n"
                          << "  Core Height (rel to origin): " << core_height << " (TR.y " << ckt_.core_region.top_right.y << " - Org.y " << origin_.y << ")\n"
                          << "  Calculated Tile W/H: (" << tile_w_ << ", " << tile_h_ << ")" << std::endl;
            }
            else
            {
                std::cout << "[Router-Warning] Cannot recalculate Tile W/H: Gx/Gy are positive, but core dimensions relative to origin are not valid for division.\n"
                          << "  Gx=" << gx_ << ", Gy=" << gy_ << ", CoreWidthRel=" << core_width << ", CoreHeightRel=" << core_height << std::endl;
                // tile_w_ and tile_h_ remain as they were (likely 0), critical error will still be hit.
            }
        }
    }
    else
    {
        std::cout << "[Router-Debug] build_grid: .route data for Gx/Gy not available or invalid (Gx="
                  << ckt_.routing_grid.grid_x << ", Gy=" << ckt_.routing_grid.grid_y
                  << "). Calculating from core region and lambda_param." << std::endl;

        origin_ = ckt_.core_region.bottom_left;
        tile_w_ = tile_h_ = lambda_param;

        if (lambda_param <= 0)
        { // Safety check if lambda_param itself is bad
            std::cout << "[Router-Warning] lambda_param (" << lambda_param << ") is non-positive. Defaulting tile size to 1.0 as a desperate measure." << std::endl;
            tile_w_ = tile_h_ = 1.0;
        }

        // Ensure tile_w_ and tile_h_ are positive before division for Gx/Gy calculation
        if (tile_w_ > 0 && tile_h_ > 0 &&
            (ckt_.core_region.top_right.x > origin_.x) &&
            (ckt_.core_region.top_right.y > origin_.y))
        {
            gx_ = (int)std::ceil((ckt_.core_region.top_right.x - origin_.x) / tile_w_);
            gy_ = (int)std::ceil((ckt_.core_region.top_right.y - origin_.y) / tile_h_);
        }
        else
        {
            gx_ = 0; // Mark Gx/Gy as invalid if tile sizes or core dimensions are not suitable
            gy_ = 0;
            std::cout << "[Router-Error] Could not determine positive tile sizes or valid core dimensions in fallback mode. Gx/Gy set to 0." << std::endl;
        }

        std::cout << "[Router-Debug] build_grid (calculated from core region):\n"
                  << "  Origin: (" << origin_.x << ", " << origin_.y << ")\n"
                  << "  Tile W/H (from lambda_param): (" << tile_w_ << ", " << tile_h_ << ") using lambda_param: " << lambda_param << "\n"
                  << "  Grid Dims (gx, gy): (" << gx_ << ", " << gy_ << ")\n"
                  << "  Core Region BL: (" << ckt_.core_region.bottom_left.x << ", " << ckt_.core_region.bottom_left.y << ")\n"
                  << "  Core Region TR: (" << ckt_.core_region.top_right.x << ", " << ckt_.core_region.top_right.y << ")" << std::endl;
    }

    // General check for invalid tile dimensions - this should be very visible
    if (tile_w_ <= 0 || tile_h_ <= 0)
    {
        std::cerr << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                  << "[Router-CRITICAL-ERROR] Invalid tile dimensions calculated in build_grid:\n"
                  << "  tile_w_ = " << tile_w_ << "\n"
                  << "  tile_h_ = " << tile_h_ << "\n"
                  << "  This will cause division by zero. Check .route file or LAMBDA constant.\n"
                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                  << std::endl;
        // It might be wise to throw an exception or std::abort() here to halt execution cleanly.
        // For example: throw std::runtime_error("Critical error: Invalid tile dimensions in Router::build_grid");
        // Or: std::abort();
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
    double div_x = (p.x - origin_.x) / tile_w_;
    double div_y = (p.y - origin_.y) / tile_h_;
    GridPt gpt = {(int)std::round(div_x),
                  (int)std::round(div_y)};

    GridPt original_gpt = gpt; // Store original for logging

    // Clamp coordinates to be within valid grid boundaries [0, gx_-1] and [0, gy_-1]
    // Ensure gx_ and gy_ are positive before trying to subtract 1 or use in std::min/max.
    // If gx_ or gy_ were somehow 0 (shouldn't be if build_grid is robust), this avoids negative upper bounds.
    if (gx_ > 0)
    {
        gpt.x = std::max(0, std::min(gpt.x, gx_ - 1));
    }
    else
    {
        gpt.x = 0; // Fallback if gx_ is not positive (e.g. 0)
    }

    if (gy_ > 0)
    {
        gpt.y = std::max(0, std::min(gpt.y, gy_ - 1));
    }
    else
    {
        gpt.y = 0; // Fallback if gy_ is not positive (e.g. 0)
    }

    // Log if clamping occurred (useful for debugging this specific scenario)
    // if (gpt.x != original_gpt.x || gpt.y != original_gpt.y)
    // {
    //     std::cout << "[Router-Debug] phys_to_grid: Clamping occurred.\n"
    //               << "  Input Point p: (" << p.x << ", " << p.y << ")\n"
    //               << "  Origin: (" << origin_.x << ", " << origin_.y << ")\n"
    //               << "  Tile W/H: (" << tile_w_ << ", " << tile_h_ << ")\n"
    //               << "  (p.x - origin_.x) / tile_w_: " << div_x << " (rounded to: " << original_gpt.x << ")\n"
    //               << "  (p.y - origin_.y) / tile_h_: " << div_y << " (rounded to: " << original_gpt.y << ")\n"
    //               << "  Original GridPt before clamp: (" << original_gpt.x << ", " << original_gpt.y << ")\n"
    //               << "  Clamped GridPt: (" << gpt.x << ", " << gpt.y << ")\n"
    //               << "  Grid Dims (gx, gy): (" << gx_ << ", " << gy_ << ")" << std::endl;
    // }

    return gpt;
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

// New method implementation
void Router::save_global_routing_results(const std::string &filename) const
{
    std::ofstream outfile(filename);
    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file for saving global routing results: " << filename << std::endl;
        return;
    }

    outfile << "Global Routing Results" << std::endl;
    outfile << "----------------------" << std::endl;

    // Save Net Paths
    outfile << std::endl
            << "Net Paths (" << net_paths_.size() << " nets routed):" << std::endl;
    for (const auto &pair : net_paths_)
    {
        const Net *net = pair.first;
        const auto &paths = pair.second;
        outfile << "Net: " << net->name << " (" << paths.size() << " path segments)" << std::endl;
        for (const auto &path_segment : paths)
        {
            outfile << "  Segment: ";
            for (size_t i = 0; i < path_segment.size(); ++i)
            {
                outfile << "(" << path_segment[i].x << "," << path_segment[i].y << ")";
                if (i < path_segment.size() - 1)
                {
                    outfile << " -> ";
                }
            }
            outfile << std::endl;
        }
    }

    // Save Horizontal Usage
    outfile << std::endl
            << "Horizontal Edge Usage (" << h_usage_.size() << " edges used):" << std::endl;
    for (const auto &pair : h_usage_)
    {
        const EdgeKey &edge = pair.first;
        int usage = pair.second;
        outfile << "Edge: [(" << edge.first.x << "," << edge.first.y << ")-(" << edge.second.x << "," << edge.second.y << ")] Usage: " << usage << std::endl;
    }

    // Save Vertical Usage
    outfile << std::endl
            << "Vertical Edge Usage (" << v_usage_.size() << " edges used):" << std::endl;
    for (const auto &pair : v_usage_)
    {
        const EdgeKey &edge = pair.first;
        int usage = pair.second;
        outfile << "Edge: [(" << edge.first.x << "," << edge.first.y << ")-(" << edge.second.x << "," << edge.second.y << ")] Usage: " << usage << std::endl;
    }

    // Save Unrouted Nets
    outfile << std::endl
            << "Unrouted Nets (" << unrouted_.size() << " nets):" << std::endl;
    for (const Net *net : unrouted_)
    {
        outfile << "Net: " << net->name << std::endl;
    }

    outfile.close();
    std::cout << "Global routing results saved to: " << filename << std::endl;
}

// Implementation for loading global routing results
std::optional<int> Router::load_global_routing_results(const std::string &filename)
{
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error: Could not open file for loading global routing results: " << filename << std::endl;
        return std::nullopt;
    }

    // Clear existing data
    net_paths_.clear();
    h_usage_.clear();
    v_usage_.clear();
    unrouted_.clear();

    std::string line, section;
    int unrouted_nets_count_from_file = -1; // To store the count read from the unrouted section header

    // Lambda to find Net* from name
    auto find_net_by_name = [&](const std::string &net_name) -> const Net *
    {
        for (const auto &pair : ckt_.net_map)
        {
            if (pair.second.name == net_name)
            {
                return &pair.second;
            }
        }
        return nullptr;
    };

    try
    {
        // Skip header lines
        std::getline(infile, line); // "Global Routing Results"
        std::getline(infile, line); // "----------------------"

        while (std::getline(infile, line))
        {
            if (line.empty())
                continue;

            if (line.find("Net Paths") != std::string::npos)
                section = "Net Paths";
            else if (line.find("Horizontal Edge Usage") != std::string::npos)
                section = "Horizontal Edge Usage";
            else if (line.find("Vertical Edge Usage") != std::string::npos)
                section = "Vertical Edge Usage";
            else if (line.find("Unrouted Nets") != std::string::npos)
            {
                section = "Unrouted Nets";
                // Parse the count from the Unrouted Nets header, e.g., "Unrouted Nets (X nets):"
                size_t open_paren = line.find('(');
                size_t space_after_count = line.find(' ', open_paren);
                if (open_paren != std::string::npos && space_after_count != std::string::npos)
                {
                    std::string count_str = line.substr(open_paren + 1, space_after_count - (open_paren + 1));
                    unrouted_nets_count_from_file = std::stoi(count_str);
                }
            }
            else if (section == "Net Paths")
            {
                if (line.rfind("Net: ", 0) == 0)
                { // Starts with "Net: "
                    std::string net_name_part = line.substr(5, line.find(" (") - 5);
                    const Net *current_net = find_net_by_name(net_name_part);
                    if (!current_net)
                    {
                        std::cerr << "Warning: Net name '" << net_name_part << "' from file not found in circuit. Skipping paths for this net." << std::endl;
                        // Continue to skip lines until next 'Net:' or new section
                        std::string temp_line;
                        std::streampos prev_pos = infile.tellg();
                        while (std::getline(infile, temp_line))
                        {
                            if (temp_line.rfind("Net: ", 0) == 0 ||
                                temp_line.find("Horizontal Edge Usage") != std::string::npos ||
                                temp_line.find("Vertical Edge Usage") != std::string::npos ||
                                temp_line.find("Unrouted Nets") != std::string::npos ||
                                temp_line.empty())
                            {
                                infile.clear();
                                infile.seekg(prev_pos);
                                break;
                            }
                            prev_pos = infile.tellg();
                        }
                        continue; // Outer loop continues
                    }

                    // Read segments for this net using peeking
                    while (infile.good())
                    {
                        std::streampos line_start_pos = infile.tellg();
                        std::string segment_line;
                        if (!std::getline(infile, segment_line))
                        {
                            break; // EOF or read error
                        }

                        if (segment_line.rfind("  Segment: ", 0) == 0)
                        {
                            std::vector<GridPt> segment;
                            std::string segment_data = segment_line.substr(11); // "(x1,y1) -> (x2,y2) ..."
                            std::istringstream segment_stream(segment_data);
                            GridPt pt;
                            char c1, comma, c2, dash, gt; // For parsing "(x,y)" and " -> "

                            // Read first point: (x,y)
                            segment_stream >> c1 >> pt.x >> comma >> pt.y >> c2;
                            if (c1 == '(' && comma == ',' && c2 == ')')
                            {
                                segment.push_back(pt);
                            }
                            else
                            {
                                std::cerr << "Warning: Malformed segment point in file: " << segment_line << std::endl;
                                continue; // Next line (hopefully a segment or new net/section)
                            }

                            // Read subsequent points: " -> (x,y)"
                            while (segment_stream >> dash >> gt >> c1 >> pt.x >> comma >> pt.y >> c2)
                            {
                                if (dash == '-' && gt == '>' && c1 == '(' && comma == ',' && c2 == ')')
                                {
                                    segment.push_back(pt);
                                }
                                else
                                {
                                    std::cerr << "Warning: Malformed segment connector or point in file: " << segment_line << std::endl;
                                    break; // Stop parsing this malformed segment
                                }
                            }
                            if (!segment.empty())
                            {
                                net_paths_[current_net].push_back(segment);
                            }
                        }
                        else
                        {
                            // Not a segment line for the current net, so put stream position back
                            infile.clear(); // Clear EOF flags if getline peeked to the end
                            infile.seekg(line_start_pos);
                            break; // Break from segment reading loop, outer loop will process this line
                        }
                    }
                    // The line that broke the inner while loop (if any) will be processed by the outer loop.
                }
            }
            else if (section == "Horizontal Edge Usage" || section == "Vertical Edge Usage")
            {
                if (line.rfind("Edge: ", 0) == 0)
                { // Starts with "Edge: "
                    GridPt p1, p2;
                    int usage;
                    // Example: Edge: [(x1,y1)-(x2,y2)] Usage: U
                    sscanf(line.c_str(), "Edge: [(%d,%d)-(%d,%d)] Usage: %d", &p1.x, &p1.y, &p2.x, &p2.y, &usage);
                    EdgeKey edge = make_edge(p1, p2);
                    if (section == "Horizontal Edge Usage")
                        h_usage_[edge] = usage;
                    else
                        v_usage_[edge] = usage;
                }
            }
            else if (section == "Unrouted Nets")
            {
                if (line.rfind("Net: ", 0) == 0)
                { // Starts with "Net: "
                    std::string net_name = line.substr(5);
                    const Net *unrouted_net = find_net_by_name(net_name);
                    if (unrouted_net)
                    {
                        unrouted_.insert(unrouted_net);
                    }
                    else
                    {
                        std::cerr << "Warning: Unrouted net name '" << net_name << "' from file not found in circuit." << std::endl;
                    }
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Exception while parsing global routing file: " << e.what() << std::endl;
        return std::nullopt;
    }

    infile.close();
    // Validate if the parsed unrouted_ set size matches the count from file if available
    if (unrouted_nets_count_from_file != -1 && (size_t)unrouted_nets_count_from_file != unrouted_.size())
    {
        std::cout << "[Router-Load] Warning: Mismatch in unrouted nets count. File header said " << unrouted_nets_count_from_file
                  << ", but loaded " << unrouted_.size() << " unrouted nets." << std::endl;
    }
    std::cout << "Successfully loaded global routing data from " << filename
              << ". Routed nets: " << net_paths_.size()
              << ", Unrouted nets: " << unrouted_.size() << std::endl;
    return unrouted_.size(); // Return the count of actually loaded unrouted nets
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

        // Defensive check before using .at()
        if (gr_.get_paths().find(&kv.second) == gr_.get_paths().end())
        {
            std::cerr << "[DetailedRouter-ERROR] Net " << kv.second.name
                      << " was not in unrouted_nets set, but also not found in get_paths() map. This indicates an issue in global router logic. Skipping net."
                      << std::endl;
            violations++; // Count this as a violation or a specific error type
            continue;
        }
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
