#ifndef VISUALIZER_H
#define VISUALIZER_H

#include "structures.h"
#include <string>
#include <unordered_map>

// Forward declaration
class Circuit;

class Visualizer
{
public:
    static void display_placement(
        const Circuit &circuit,
        const std::unordered_map<std::string, Point> &placement,
        const std::string &output_file = "placement.svg", // Default to SVG
        bool show_labels = false                          // Added flag for labels
    );

private:
    // Prevent instantiation of this utility class
    Visualizer() = delete;

    // Helper function for SVG output
    static void generate_svg(const Circuit &circuit,
                             const std::unordered_map<std::string, Point> &placement,
                             const std::string &filename,
                             bool show_labels = false); // Added flag for labels
};

#endif // VISUALIZER_H
