#include "../include/parser.h"
#include "../include/placer.h"
#include "../include/fm.h"
#include "../include/visualizer.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iomanip>

void print_usage(const char *program_name)
{
    std::cout << "Usage: " << program_name << " [options] <dataset_name>\n"
              << "  <dataset_name>: Name of the benchmark dataset (e.g., superblue18). Assumes files are in ./data/<dataset_name>/<dataset_name>.*\n"
              << "Options:\n"
              << "  -h, --help              Show this help message\n"
              << "  -s, --strategy <type>   Placement strategy (default: bisection)\n"
              << "                          Available strategies:\n"
              << "                            - bisection\n"
              << "                            - quadrature\n"
              << "                            - slice-bisection\n"
              << "                            - cut-oriented\n"
              << "  -o, --output <file>     Output file for placement results (default: placement.out)\n"
              << "  --svg <file>          Output file for SVG visualization (default: final_placement.svg)\n"
              << "  -v, --verbose           Enable verbose output\n"
              << "\n" // Add a newline for spacing
              << "Example:\n"
              << "  " << program_name << " -s quadrature -o my_placement.txt --svg my_layout.svg superblue18\n"
              << std::endl;
}

PlacementStrategy parse_strategy(const std::string &strategy_str)
{
    if (strategy_str == "bisection")
        return BISECTION;
    if (strategy_str == "quadrature")
        return QUADRATURE;
    if (strategy_str == "slice-bisection")
        return SLICE_BISECTION;
    if (strategy_str == "cut-oriented")
        return CUT_ORIENTED;

    throw std::runtime_error("Unknown placement strategy: " + strategy_str);
}

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Default values
    std::string dataset_name;
    std::string strategy_str = "bisection";
    std::string output_file = "placement.out";
    std::string svg_output_file = "final_placement.svg";
    bool verbose = false;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help")
        {
            print_usage(argv[0]);
            return 0;
        }
        else if (arg == "-s" || arg == "--strategy")
        {
            if (++i < argc)
            {
                strategy_str = argv[i];
            }
            else
            {
                std::cerr << "Error: Strategy argument requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        else if (arg == "-o" || arg == "--output")
        {
            if (++i < argc)
            {
                output_file = argv[i];
            }
            else
            {
                std::cerr << "Error: Output argument requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--svg")
        {
            if (++i < argc)
            {
                svg_output_file = argv[i];
            }
            else
            {
                std::cerr << "Error: SVG output argument requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        else if (arg == "-v" || arg == "--verbose")
        {
            verbose = true;
        }
        else
        {
            // Assume the last argument is the dataset name
            if (!dataset_name.empty())
            {
                std::cerr << "Error: Multiple dataset names provided or dataset name specified before options." << std::endl;
                print_usage(argv[0]);
                return 1;
            }
            dataset_name = arg;
        }
    }

    // Check if dataset name is provided
    if (dataset_name.empty())
    {
        std::cerr << "Error: Benchmark dataset name is required" << std::endl;
        print_usage(argv[0]);
        return 1;
    }

    // Construct paths
    std::string benchmark_dir = "./data/" + dataset_name;
    std::string aux_file_path = benchmark_dir + "/" + dataset_name + ".aux";

    try
    {
        // Create parser and circuit
        Parser parser;
        Circuit circuit;
        // Region placable_region;

        // Parse benchmark files
        std::cout << "Parsing benchmark files for dataset: " << dataset_name << " from directory: " << benchmark_dir << std::endl;
        parser.parse_ispd2011(aux_file_path, circuit);

        // --- Generate Initial Visualization ---
        std::cout << "Generating initial placement visualization..." << std::endl;
        std::map<std::string, Point> initial_placement;
        for (const auto &[name, node] : circuit.cell_map)
        {
            initial_placement[name] = node.pos; // Use initial positions from node data
        }
        try
        {
            Visualizer::display_placement(circuit, initial_placement, "initial_placement.svg");
            std::cout << "Initial placement visualization saved to: initial_placement.svg" << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Warning: Could not generate initial SVG visualization: " << e.what() << std::endl;
        }
        // --------------------------------------

        // Create FM partitioner
        FMPartitioner fm(circuit);

        // Create placer
        Placer placer(circuit, fm);

        // Parse placement strategy
        PlacementStrategy strategy = parse_strategy(strategy_str);

        // Run placement
        std::cout << "Starting placement with strategy: " << strategy_str << std::endl;
        placer.place(strategy);

        // Get placement results
        const auto &placement = placer.get_placement();

        // --- Generate Final Visualization ---
        std::cout << "Generating final placement visualization..." << std::endl;
        try
        {
            Visualizer::display_placement(circuit, placement, svg_output_file);
            std::cout << "Final placement visualization saved to: " << svg_output_file << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Warning: Could not generate final SVG visualization: " << e.what() << std::endl;
        }
        // ------------------------------------

        // Write results to output file
        std::ofstream out(output_file);
        if (!out)
        {
            throw std::runtime_error("Failed to open output file: " + output_file);
        }

        // Write placement results in ISPD format
        for (const auto &[node_name, pos] : placement)
        {
            out << node_name << " " << pos.x << " " << pos.y << std::endl;
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Total execution time: " << std::fixed << std::setprecision(3) << duration.count() / 1000.0 << " seconds." << std::endl;
        return 1;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Placement completed successfully. Results written to: " << output_file << std::endl;
    std::cout << "Total execution time: " << std::fixed << std::setprecision(3) << duration.count() / 1000.0 << " seconds." << std::endl;

    return 0;
}