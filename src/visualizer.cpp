#include "../include/visualizer.h"
#include "../include/structures.h" // Contains Circuit, Node, Point, Rect definitions
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <limits>    // For numeric_limits
#include <algorithm> // For std::min/max
#include <cmath>     // For std::fabs

// Helper function to generate SVG output
void Visualizer::generate_svg(const Circuit &circuit,
                              const std::map<std::string, Point> &placement,
                              const std::string &filename)
{
    std::ofstream svg_file(filename);
    if (!svg_file)
    {
        throw std::runtime_error("Error: Cannot open SVG output file: " + filename);
    }

    // --- Calculate Bounding Box ---
    double min_x = circuit.core_region.bottom_left.x;
    double min_y = circuit.core_region.bottom_left.y;
    double max_x = circuit.core_region.top_right.x;
    double max_y = circuit.core_region.top_right.y + circuit.core_rows[circuit.core_rows.size() - 1].height;

    if (placement.empty() && (max_x == min_x || max_y == min_y))
    {
        // Handle empty placement or zero-area core region
        min_x = 0;
        min_y = 0;
        max_x = 100;
        max_y = 100; // Default view
    }
    else if (!placement.empty())
    {
        // Expand bounds based on placed nodes
        for (const auto &[node_name, pos] : placement)
        {
            try
            {
                const Node &node = circuit.cell_map.at(node_name); // Use .at() for bounds checking
                min_x = std::min(min_x, pos.x);
                min_y = std::min(min_y, pos.y);
                max_x = std::max(max_x, pos.x + node.width);
                max_y = std::max(max_y, pos.y + node.height);
            }
            catch (const std::out_of_range &oor)
            {
                // This should ideally not happen if placement map is consistent with circuit.cell_map
                // Consider logging a warning here if necessary.
                continue;
            }
        }
    }

    // Add padding
    double range_x = max_x - min_x;
    double range_y = max_y - min_y;
    double padding = 0.05 * std::max(range_x > 1e-9 ? range_x : 100.0, range_y > 1e-9 ? range_y : 100.0); // 5% padding, handle zero range
    if (padding < 1e-9)
        padding = 10.0; // Minimum padding if range is very small

    double view_min_x = min_x - padding;
    double view_min_y = min_y - padding;
    double view_width = (max_x + padding) - view_min_x;
    double view_height = (max_y + padding) - view_min_y;

    // Ensure non-zero dimensions for viewBox
    if (view_width < 1e-9)
        view_width = 100.0;
    if (view_height < 1e-9)
        view_height = 100.0;

    // --- Write SVG Header ---
    // Fixed SVG size, viewBox handles scaling and aspect ratio
    svg_file << R"(<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="1000")"
             << R"( viewBox=")" << view_min_x << " " << view_min_y << " "
             << view_width << " " << view_height << R"(" preserveAspectRatio="xMidYMid meet">)" << std::endl;

    // --- Draw SVG Border ---
    svg_file << R"(  <rect x=")" << view_min_x << R"(" y=")" << view_min_y
             << R"(" width=")" << view_width << R"(" height=")" << view_height
             << R"(" fill="none" stroke="green" stroke-width="1" />)" << std::endl; // Added 1px border

    // --- Add Coordinate System Transformation ---
    // Flips Y-axis and translates origin to match standard Cartesian coordinates (origin bottom-left)
    // Using R"svg(...)svg" to avoid issues with parentheses inside the transform attribute
    svg_file << R"svg(  <g transform="translate(0, )svg" << (view_min_y + view_height) << R"svg() scale(1, -1)">)svg" << std::endl;

    // Determine a reasonable stroke width based on viewbox size
    double stroke_width = std::max(view_width, view_height) / 1000.0;
    if (stroke_width < 1e-9)
        stroke_width = 0.1; // Minimum stroke width

    // --- Draw Core Region ---
    double core_width = circuit.core_region.top_right.x - circuit.core_region.bottom_left.x;
    double core_height = circuit.core_region.top_right.y - circuit.core_region.bottom_left.y;
    svg_file << R"(    <rect x=")" << circuit.core_region.bottom_left.x
             << R"(" y=")" << circuit.core_region.bottom_left.y
             << R"(" width=")" << core_width
             << R"(" height=")" << core_height
             << R"(" fill="none" stroke="red" stroke-width=")" << stroke_width << R"(" />)" << std::endl;

    // --- Draw Placed Nodes ---
    for (const auto &[node_name, pos] : placement)
    {
        try
        {
            const Node &node = circuit.cell_map.at(node_name); // Use .at() for safety
            std::string fill_color = "lightblue";              // Default for movable
            if (node.type == FIXED || node.type == TERMINAL || node.type == TERMINAL_NI)
            {
                fill_color = "lightcoral"; // Fixed/Terminal nodes
            }

            svg_file << R"(    <rect x=")" << pos.x
                     << R"(" y=")" << pos.y
                     << R"(" width=")" << node.width
                     << R"(" height=")" << node.height
                     << R"(" fill=")" << fill_color
                     << R"(" stroke="blue" stroke-width=")" << stroke_width / 2.0 << R"(" />)" << std::endl;

            // Optional: Add text label (can clutter the view for many nodes)
            // double text_x = pos.x + node.width / 2.0;
            // double text_y = pos.y + node.height / 2.0;
            // double font_size = std::min(node.width, node.height) / 5.0; // Example sizing
            // if (font_size > 1e-9) {
            //      svg_file << R"(    <text x=")" << text_x << R"(" y=")" << text_y
            //               << R"(" font-size=")" << font_size << R"(" text-anchor="middle" dominant-baseline="middle")"
            //                // Apply inverse transform to text so it's readable
            //               << R"( transform="translate(0, )" << 2 * text_y << R"() scale(1, -1))">)"
            //               << node.name << R"(</text>)" << "
            // "; // This semicolon was causing an error
            // }
        }
        catch (const std::out_of_range &oor)
        {
            // Node from placement map not found in cell_map, skip drawing
            continue;
        }
    }

    // --- Close Transformation Group ---
    svg_file << "  </g>" << std::endl;

    // --- Write SVG Footer ---
    svg_file << "</svg>" << std::endl;

    svg_file.close();
}

// Public method implementation
void Visualizer::display_placement(const Circuit &circuit,
                                   const std::map<std::string, Point> &placement,
                                   const std::string &output_file)
{
    // Currently, just calls the SVG generation helper.
    // Could be extended to choose different output formats later.
    generate_svg(circuit, placement, output_file);
}