#ifndef VISUALIZER_H
#define VISUALIZER_H

#include "structures.h"
#include <string>
#include <map>

// Forward declaration
class Circuit;

class Visualizer
{
public:
    // Static method to generate a visual representation of the placement.
    // Could output an SVG, PNG (using a library like Cairo or SDL),
    // or data for an external tool like gnuplot or Python/matplotlib.
    static void display_placement(
        const Circuit &circuit,
        const std::map<std::string, Point> &placement,
        const std::string &output_file = "placement.svg" // Default to SVG
    );

    // Optional: Add methods to visualize congestion maps, specific nets, etc.
    // static void display_congestion_map(...);
    // static void display_net(...);

private:
    // Prevent instantiation of this utility class
    Visualizer() = delete;

    // Helper function for SVG output (example)
    static void generate_svg(const Circuit &circuit,
                             const std::map<std::string, Point> &placement,
                             const std::string &filename);

    // Add helpers for other output formats if needed
    // static void generate_png(...);
    // static void generate_gnuplot_data(...);
};

#endif // VISUALIZER_H
