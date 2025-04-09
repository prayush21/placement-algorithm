#include "../include/parser.h"
#include "../include/placer.h"
#include "../include/fm.h"
#include "../include/evaluator.h"
#include "../include/visualizer.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <unordered_map>
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
              << "  -a, --all               Run all strategies\n"
              << "  -o, --output <file>     Output file for placement results (default: <dataset>_<strategy>_pl.txt)\n"
              << "  --svg <file>          Output file for SVG visualization (default: <dataset>_<strategy>_final.svg)\n"
              << "  --eval <file>         Output file for evaluation metrics (default: <dataset>_<strategy>_eval.txt)\n"
              << "\n" // Add a newline for spacing
              << "Examples:\n"
              << "  " << program_name << " -s quadrature -o my_placement.txt --svg my_layout.svg --eval my_eval.txt superblue18\n"
              << "  " << program_name << " -s bisection superblue18\n"
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

void run_multiple_strategies(const std::string &dataset_name)
{
    std::vector<std::string> strategies = {"quadrature", "bisection", "slice-bisection"};

    std::cout << "Starting placement runs with different strategies..." << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (const auto &strategy : strategies)
    {
        std::cout << "Running placement with strategy: " << strategy << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        // Create new output filenames for each strategy
        std::string output_file = dataset_name + "_" + strategy + "_pl.txt";
        std::string svg_output_file = dataset_name + "_" + strategy + "_final.svg";
        std::string evaluation_output_file = dataset_name + "_" + strategy + "_eval.txt";

        // Construct paths
        std::string benchmark_dir = "./data/" + dataset_name;
        std::string aux_file_path = benchmark_dir + "/" + dataset_name + ".aux";

        try
        {
            // Create parser and circuit
            Parser parser;
            Circuit circuit;

            // Parse benchmark files
            std::cout << "Parsing benchmark files for dataset: " << dataset_name << " from directory: " << benchmark_dir << std::endl;
            parser.parse_ispd2011(aux_file_path, circuit);

            // Generate initial visualization
            std::cout << "Generating initial placement visualization..." << std::endl;
            std::unordered_map<std::string, Point> initial_placement;
            for (const auto &[name, node] : circuit.cell_map)
            {
                initial_placement[name] = node.pos;
            }
            try
            {
                std::string initial_svg_file = dataset_name + "_" + strategy + "_initial.svg";
                Visualizer::display_placement(circuit, initial_placement, initial_svg_file, false);
                std::cout << "Initial placement visualization saved to: " << initial_svg_file << std::endl;
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Could not generate initial SVG visualization: " << e.what() << std::endl;
            }

            // Create FM partitioner and placer
            FMPartitioner fm(circuit);
            Placer placer(circuit, fm);

            // Parse placement strategy
            PlacementStrategy placement_strategy = parse_strategy(strategy);

            // Run placement
            std::cout << "Starting placement with strategy: " << strategy << std::endl;
            placer.place(placement_strategy);

            // Get placement results
            const auto &placement = placer.get_placement();

            // Generate final visualization
            std::cout << "Generating final placement visualization..." << std::endl;
            try
            {
                Visualizer::display_placement(circuit, placement, svg_output_file, false);
                std::cout << "Final placement visualization saved to: " << svg_output_file << std::endl;
            }
            catch (const std::exception &e)
            {
                std::cerr << "Warning: Could not generate final SVG visualization: " << e.what() << std::endl;
            }

            // Write results to output file
            std::ofstream out(output_file);
            if (!out)
            {
                throw std::runtime_error("Failed to open output file: " + output_file);
            }

            // Write placement results
            for (const auto &[node_name, pos] : placement)
            {
                out << node_name << " " << pos.x << " " << pos.y << std::endl;
            }
            out.close();

            // Calculate and write evaluation metrics
            std::cout << "Calculating evaluation metrics..." << std::endl;
            Evaluator evaluator(circuit, placement);
            double total_hpwl = evaluator.calculate_total_hpwl();

            double total_movable_area = 0.0;
            for (const auto &[node_name, pos] : placement)
            {
                auto it = circuit.cell_map.find(node_name);
                if (it != circuit.cell_map.end() && it->second.type == MOVABLE)
                {
                    total_movable_area += it->second.area;
                }
            }

            std::ofstream eval_out(evaluation_output_file);
            if (!eval_out)
            {
                throw std::runtime_error("Failed to open evaluation output file: " + evaluation_output_file);
            }
            eval_out << "Dataset: " << dataset_name << std::endl;
            eval_out << "Placement Strategy: " << strategy << std::endl;
            eval_out << "Total Movable Node Area: " << std::fixed << std::setprecision(2) << total_movable_area << std::endl;
            eval_out << "Total HPWL: " << std::fixed << std::setprecision(2) << total_hpwl << std::endl;
            eval_out.close();

            std::cout << "----------------------------------------" << std::endl;
            std::cout << "Completed placement with strategy: " << strategy << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error running strategy " << strategy << ": " << e.what() << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << "Failed placement with strategy: " << strategy << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            std::cout << std::endl;
        }
    }

    std::cout << "All placement runs completed!" << std::endl;
}

int main(int argc, char *argv[])
{
    auto start_time = std::chrono::high_resolution_clock::now();

    // Default values
    std::string dataset_name;
    std::string strategy_str = "bisection";
    std::string output_file;
    std::string svg_output_file;
    std::string evaluation_output_file;
    bool run_all_strategies = false;

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
        else if (arg == "-a" || arg == "--all")
        {
            run_all_strategies = true;
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
        else if (arg == "--eval")
        {
            if (++i < argc)
            {
                evaluation_output_file = argv[i];
            }
            else
            {
                std::cerr << "Error: Evaluation output argument requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
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

    if (run_all_strategies)
    {
        run_multiple_strategies(dataset_name);
        return 0;
    }

    // Set default filenames if not provided
    if (output_file.empty())
    {
        output_file = dataset_name + "_" + strategy_str + "_pl.txt";
    }
    if (svg_output_file.empty())
    {
        svg_output_file = dataset_name + "_" + strategy_str + "_final.svg";
    }
    if (evaluation_output_file.empty())
    {
        evaluation_output_file = dataset_name + "_" + strategy_str + "_eval.txt";
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
        std::unordered_map<std::string, Point> initial_placement;
        for (const auto &[name, node] : circuit.cell_map)
        {
            initial_placement[name] = node.pos; // Use initial positions from node data
        }
        try
        {
            std::string initial_svg_file = dataset_name + "_" + strategy_str + "_initial.svg";
            Visualizer::display_placement(circuit, initial_placement, initial_svg_file, false);
            std::cout << "Initial placement visualization saved to: " << initial_svg_file << std::endl;
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
            Visualizer::display_placement(circuit, placement, svg_output_file, false);
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
        out.close();

        // --- Calculate Evaluation Metrics ---
        std::cout << "Calculating evaluation metrics..." << std::endl;
        Evaluator evaluator(circuit, placement);
        double total_hpwl = evaluator.calculate_total_hpwl();

        double total_movable_area = 0.0;
        for (const auto &[node_name, pos] : placement)
        {
            auto it = circuit.cell_map.find(node_name);
            if (it != circuit.cell_map.end() && it->second.type == MOVABLE)
            {
                total_movable_area += it->second.area;
            }
        }

        auto end_time_eval = std::chrono::high_resolution_clock::now();
        auto duration_eval = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_eval - start_time);
        double execution_time_sec = duration_eval.count() / 1000.0;

        // --- Write Evaluation File ---
        std::cout << "Writing evaluation metrics to: " << evaluation_output_file << std::endl;
        std::ofstream eval_out(evaluation_output_file);
        if (!eval_out)
        {
            throw std::runtime_error("Failed to open evaluation output file: " + evaluation_output_file);
        }
        eval_out << "Dataset: " << dataset_name << std::endl;
        eval_out << "Placement Strategy: " << strategy_str << std::endl;
        eval_out << "Execution Time (seconds): " << std::fixed << std::setprecision(3) << execution_time_sec << std::endl;
        eval_out << "Total Movable Node Area: " << std::fixed << std::setprecision(2) << total_movable_area << std::endl;
        eval_out << "Total HPWL: " << std::fixed << std::setprecision(2) << total_hpwl << std::endl;
        eval_out.close();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        auto end_time_err = std::chrono::high_resolution_clock::now();
        auto duration_err = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_err - start_time);
        std::cout << "Total execution time (before error): " << std::fixed << std::setprecision(3) << duration_err.count() / 1000.0 << " seconds." << std::endl;
        return 1;
    }

    return 0;
}