#include "../include/parser.h"
#include "../include/placer.h"
#include "../include/fm.h"
#include "../include/router.h"
#include "../include/evaluator.h"
#include "../include/visualizer.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>

// Routing algorithm enums
// enum GlobalAlgo
// {
//     ALGO_HADLOCK,
//     ALGO_SOUKUP
// };

// enum MergeAlgo
// {
//     MERGE_STAR,
//     MERGE_MST
// };

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
              << "  --router <type>       Global routing algorithm (default: hadlock)\n"
              << "                          Available algorithms:\n"
              << "                            - hadlock\n"
              << "                            - soukup\n"
              << "  --route_from <file>   Path to an existing placement file (.pl.txt format). If provided, skips placement and routes from this file. Cannot be used with -s or -a.\n"
              << "\n" // Add a newline for spacing
              << "Examples:\n"
              << "  " << program_name << " -s quadrature -o my_placement.txt --svg my_layout.svg --eval my_eval.txt superblue18\n"
              << "  " << program_name << " -s bisection superblue18\n"
              << "  " << program_name << " --route_from existing_placement.pl.txt --router soukup superblue18\n"
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

// Helper function to load placement from a file
bool load_placement_from_file(const std::string &filepath, Circuit &circuit, std::unordered_map<std::string, Point> &placement_map)
{
    std::ifstream in_file(filepath);
    if (!in_file)
    {
        std::cerr << "Error: Cannot open placement file: " << filepath << std::endl;
        return false;
    }
    std::string line;
    int line_num = 0;
    while (std::getline(in_file, line))
    {
        line_num++;
        std::istringstream iss(line);
        std::string cell_name;
        double x, y;

        // Skip empty lines or comment lines (often starting with # or specific keywords)
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        // A simple check for header-like lines if they contain typical delimiters and no numbers at start
        if (line.find(":") != std::string::npos && !(std::isdigit(line[0]) || line[0] == '-'))
        {
            std::cout << "Skipping potential header line in placement file: " << line << std::endl;
            continue;
        }

        if (!(iss >> cell_name >> x >> y))
        {
            std::cerr << "Warning: Malformed line #" << line_num << " in placement file: \"" << line << "\". Expected format: cell_name x y" << std::endl;
            continue; // Skip malformed lines
        }

        auto it = circuit.cell_map.find(cell_name);
        if (it != circuit.cell_map.end())
        {
            it->second.pos = {x, y};
            placement_map[cell_name] = {x, y};
        }
        else
        {
            std::cerr << "Warning: Cell '" << cell_name << "' from placement file (line " << line_num << ") not found in circuit definition. Skipping." << std::endl;
        }
    }
    if (placement_map.empty() && line_num > 0)
    {
        std::cerr << "Warning: No valid cell placements were loaded from file: " << filepath << ". Ensure the file is not empty and cells match the dataset." << std::endl;
    }
    else if (placement_map.empty())
    {
        std::cerr << "Warning: Placement file " << filepath << " seems to be empty or contains no valid data." << std::endl;
    }
    return true;
}

// Function to run a single placement/routing flow
void run_single_flow(
    const std::string &dataset_name,
    const std::string &strategy_str_param,
    const std::string &router_algo_str,
    const std::string &placement_output_file_param,
    const std::string &svg_output_file_param,
    const std::string &evaluation_output_file_param,
    const std::string &route_from_filepath_param,
    std::chrono::time_point<std::chrono::high_resolution_clock> flow_start_time)
{
    std::string current_strategy_label = strategy_str_param;

    // Construct paths
    std::string benchmark_dir = "./data/" + dataset_name;
    std::string aux_file_path = benchmark_dir + "/" + dataset_name + ".aux";

    Circuit circuit;
    Parser parser;

    std::cout << "Parsing benchmark files for dataset: " << dataset_name << " from directory: " << benchmark_dir << std::endl;
    parser.parse_ispd2011(aux_file_path, circuit);

    std::unordered_map<std::string, Point> placement_to_use;

    if (!route_from_filepath_param.empty())
    {
        std::cout << "Loading placement from file: " << route_from_filepath_param << std::endl;
        if (!load_placement_from_file(route_from_filepath_param, circuit, placement_to_use))
        {
            throw std::runtime_error("Failed to load placement from file: " + route_from_filepath_param);
        }
        // Use a descriptive label for "strategy" when loading from file
        current_strategy_label = "custom_from_" + std::filesystem::path(route_from_filepath_param).stem().string();

        std::cout << "Generating visualization of loaded placement..." << std::endl;
        try
        {
            // Ensure unique filename for this specific visualization
            std::string loaded_svg_file = dataset_name + "_" + current_strategy_label + "_" + router_algo_str + "_loaded_placement.svg";
            Visualizer::display_placement(circuit, placement_to_use, loaded_svg_file, false);
            std::cout << "Loaded placement visualization saved to: " << loaded_svg_file << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Warning: Could not generate loaded placement SVG: " << e.what() << std::endl;
        }
    }
    else
    {
        std::cout << "Generating initial placement visualization (from parser)..." << std::endl;
        std::unordered_map<std::string, Point> initial_placement_from_parser;
        for (const auto &[name, node] : circuit.cell_map)
        {
            initial_placement_from_parser[name] = node.pos;
        }
        try
        {
            std::string initial_svg_file = dataset_name + "_" + current_strategy_label + "_" + router_algo_str + "_initial_parsed.svg";
            Visualizer::display_placement(circuit, initial_placement_from_parser, initial_svg_file, true);
            std::cout << "Initial (parsed) placement visualization saved to: " << initial_svg_file << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Warning: Could not generate initial (parsed) SVG: " << e.what() << std::endl;
        }

        FMPartitioner fm(circuit);
        Placer placer(circuit, fm);
        PlacementStrategy placement_strategy_enum = parse_strategy(strategy_str_param); // Use original strategy_str_param here

        std::cout << "Starting placement with strategy: " << strategy_str_param << std::endl;
        placer.place(placement_strategy_enum);
        placement_to_use = placer.get_placement();

        if (!placement_output_file_param.empty())
        {
            std::cout << "Writing placement results to: " << placement_output_file_param << std::endl;
            std::ofstream out(placement_output_file_param);
            if (!out)
            {
                throw std::runtime_error("Failed to open placement output file: " + placement_output_file_param);
            }
            for (const auto &[node_name, pos] : placement_to_use)
            {
                out << node_name << " " << pos.x << " " << pos.y << std::endl;
            }
            out.close();
        }
        else
        {
            std::cout << "Note: No placement output file specified (-o). Results will not be saved to a .pl.txt file." << std::endl;
        }
    }

    // Routing and Evaluation common part
    Router grouter(circuit);
    GlobalAlgo router_algo_enum = (router_algo_str == "hadlock") ? ALGO_HADLOCK : ALGO_SOUKUP;
    std::cout << "Starting global routing with algorithm: " << router_algo_str << std::endl;
    int unrouted_nets = grouter.global_route(router_algo_enum, false, MERGE_STAR);

    DetailedRouter dr(circuit, grouter);
    std::cout << "Starting detailed routing..." << std::endl;
    int drc_violations = dr.detailed_route();

    std::cout << "Generating final layout visualization..." << std::endl;
    try
    {
        Visualizer::display_placement(circuit, placement_to_use, svg_output_file_param, false);
        std::cout << "Final layout visualization saved to: " << svg_output_file_param << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Warning: Could not generate final SVG: " << e.what() << std::endl;
    }

    std::cout << "Calculating evaluation metrics..." << std::endl;
    Evaluator evaluator(circuit, placement_to_use);
    double total_hpwl = evaluator.calculate_total_hpwl();

    double total_movable_area = 0.0;
    for (const auto &[node_name, node_data] : circuit.cell_map)
    {
        if (node_data.type == MOVABLE)
        {
            total_movable_area += node_data.area;
        }
    }

    auto end_time_eval = std::chrono::high_resolution_clock::now();
    auto duration_eval = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_eval - flow_start_time);
    double execution_time_sec = duration_eval.count() / 1000.0;

    std::ofstream eval_out(evaluation_output_file_param);
    if (!eval_out)
    {
        throw std::runtime_error("Failed to open evaluation output file: " + evaluation_output_file_param);
    }
    eval_out << "Dataset: " << dataset_name << std::endl;
    if (!route_from_filepath_param.empty())
    {
        eval_out << "Mode: Routed from existing placement file" << std::endl;
        eval_out << "Input Placement File: " << route_from_filepath_param << std::endl;
    }
    else
    {
        eval_out << "Mode: Full Place and Route" << std::endl;
        eval_out << "Placement Strategy: " << strategy_str_param << std::endl;
    }
    eval_out << "Routing Algorithm: " << router_algo_str << std::endl;
    eval_out << "Execution Time (seconds): " << std::fixed << std::setprecision(3) << execution_time_sec << std::endl;
    eval_out << "Total Movable Node Area: " << std::fixed << std::setprecision(2) << total_movable_area << std::endl;
    eval_out << "Total HPWL: " << std::fixed << std::setprecision(2) << total_hpwl << std::endl;
    eval_out << "Unrouted Nets (Global): " << unrouted_nets << std::endl;
    eval_out << "DRC Violations (Detailed): " << drc_violations << std::endl;
    eval_out.close();

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Completed flow. ";
    if (!route_from_filepath_param.empty())
    {
        std::cout << "Routed from: " << route_from_filepath_param;
    }
    else
    {
        std::cout << "Strategy: " << strategy_str_param;
    }
    std::cout << ", Router: " << router_algo_str << std::endl;
    std::cout << "Evaluation results saved to: " << evaluation_output_file_param << std::endl;
    std::cout << "----------------------------------------" << std::endl;
}

void run_multiple_strategies(const std::string &dataset_name, const std::string &cli_router_algo_str)
{
    std::vector<std::string> strategies = {"quadrature", "bisection", "slice-bisection", "cut-oriented"};
    // Use router algorithm from CLI if provided, otherwise default to hadlock for --all runs if not specified.
    // The main function will pass its router_algo_str here.

    std::cout << "Starting placement runs with different strategies using router: " << cli_router_algo_str << "..." << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    for (const auto &strategy : strategies)
    {
        std::cout << "Running placement with strategy: " << strategy << std::endl;
        std::cout << "----------------------------------------" << std::endl;

        // Create new output filenames for each strategy
        std::string output_file = dataset_name + "_" + strategy + "_" + cli_router_algo_str + "_pl.txt";
        std::string svg_output_file = dataset_name + "_" + strategy + "_" + cli_router_algo_str + "_final.svg";
        std::string evaluation_output_file = dataset_name + "_" + strategy + "_" + cli_router_algo_str + "_eval.txt";

        // Construct paths
        std::string benchmark_dir = "./data/" + dataset_name;
        std::string aux_file_path = benchmark_dir + "/" + dataset_name + ".aux";

        try
        {
            auto start_time = std::chrono::high_resolution_clock::now();
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
                std::string initial_svg_file = dataset_name + "_" + strategy + "_" + cli_router_algo_str + "_initial.svg";
                Visualizer::display_placement(circuit, initial_placement, initial_svg_file, true);
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

            // Create router and run global routing
            Router grouter(circuit);
            GlobalAlgo router_algo = (cli_router_algo_str == "hadlock") ? ALGO_HADLOCK : ALGO_SOUKUP;
            int unrouted_nets = grouter.global_route(router_algo, false, MERGE_STAR);

            DetailedRouter dr(circuit, grouter);
            int drc_violations = dr.detailed_route(); // returns overflow count etc.

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

            auto end_time_eval = std::chrono::high_resolution_clock::now();
            auto duration_eval = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_eval - start_time);
            double execution_time_sec = duration_eval.count() / 1000.0;

            // Write evaluation metrics to file
            std::ofstream eval_out(evaluation_output_file);
            if (!eval_out)
            {
                throw std::runtime_error("Failed to open evaluation output file: " + evaluation_output_file);
            }
            eval_out << "Dataset: " << dataset_name << std::endl;
            eval_out << "Placement Strategy: " << strategy << std::endl;
            eval_out << "Execution Time (seconds): " << std::fixed << std::setprecision(3) << execution_time_sec << std::endl;
            eval_out << "Total Movable Node Area: " << std::fixed << std::setprecision(2) << total_movable_area << std::endl;
            eval_out << "Total HPWL: " << std::fixed << std::setprecision(2) << total_hpwl << std::endl;
            eval_out << "Unrouted Nets: " << unrouted_nets << std::endl;
            eval_out << "DRC Violations: " << drc_violations << std::endl;
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
    auto program_start_time = std::chrono::high_resolution_clock::now();

    // Default values
    std::string dataset_name;
    std::string strategy_str = "bisection";  // Default if -s not used
    std::string router_algo_str = "hadlock"; // Default router algorithm
    bool run_all_strategies_flag = false;
    std::string cli_placement_output_file;  // From -o
    std::string cli_svg_output_file;        // From --svg
    std::string cli_evaluation_output_file; // From --eval
    std::string cli_route_from_filepath;    // From --route_from

    // Parse command line arguments
    bool strategy_explicitly_set = false;
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
                strategy_explicitly_set = true;
            }
            else
            {
                std::cerr << "Error: Strategy argument requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--router")
        {
            if (++i < argc)
            {
                router_algo_str = argv[i];
            }
            else
            {
                std::cerr << "Error: Router algorithm not specified" << std::endl;
                return 1;
            }
            // Validate router algorithm
            if (router_algo_str != "hadlock" && router_algo_str != "soukup")
            {
                std::cerr << "Error: Invalid router algorithm. Must be 'hadlock' or 'soukup'" << std::endl;
                return 1;
            }
        }
        else if (arg == "-a" || arg == "--all")
        {
            run_all_strategies_flag = true;
        }
        else if (arg == "-o" || arg == "--output")
        {
            if (++i < argc)
            {
                cli_placement_output_file = argv[i];
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
                cli_svg_output_file = argv[i];
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
                cli_evaluation_output_file = argv[i];
            }
            else
            {
                std::cerr << "Error: Evaluation output argument requires a value" << std::endl;
                print_usage(argv[0]);
                return 1;
            }
        }
        else if (arg == "--route_from")
        {
            if (++i < argc)
            {
                cli_route_from_filepath = argv[i];
            }
            else
            {
                std::cerr << "Error: --route_from argument requires a filepath value" << std::endl;
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

    // Validations for --route_from
    if (!cli_route_from_filepath.empty())
    {
        if (run_all_strategies_flag)
        {
            std::cerr << "Error: --route_from cannot be used with -a/--all." << std::endl;
            print_usage(argv[0]);
            return 1;
        }
        if (strategy_explicitly_set)
        {
            std::cerr << "Error: --route_from cannot be used with -s/--strategy." << std::endl;
            print_usage(argv[0]);
            return 1;
        }
        if (!cli_placement_output_file.empty())
        {
            std::cout << "Warning: -o/--output (for placement file) is ignored when --route_from is used, as placement is not generated." << std::endl;
            // No need to clear cli_placement_output_file, run_single_flow won't use it if cli_route_from_filepath is set
        }
    }

    if (run_all_strategies_flag)
    {
        run_multiple_strategies(dataset_name, router_algo_str); // Pass router_algo_str
        // Calculate total time for --all if needed, or rely on individual timings
        auto total_end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(total_end_time - program_start_time);
        std::cout << "Total execution time for all strategies: " << std::fixed << std::setprecision(3) << total_duration.count() / 1000.0 << " seconds." << std::endl;
        return 0;
    }

    // Determine effective names for output files for a single run
    std::string effective_label_for_output;
    if (!cli_route_from_filepath.empty())
    {
        std::string input_stem = std::filesystem::path(cli_route_from_filepath).stem().string();
        // To avoid excessively long filenames, perhaps a simpler label
        effective_label_for_output = "from_" + input_stem;
        // Or just "custom_placement"
        // effective_label_for_output = "custom_routed";
    }
    else
    {
        effective_label_for_output = strategy_str;
    }

    std::string final_placement_output_file = cli_placement_output_file;
    if (final_placement_output_file.empty() && cli_route_from_filepath.empty())
    { // Only default if not routing from file and not specified
        final_placement_output_file = dataset_name + "_" + effective_label_for_output + "_" + router_algo_str + "_pl.txt";
    }
    else if (final_placement_output_file.empty() && !cli_route_from_filepath.empty())
    {
        // No default placement output file if --route_from is used
        final_placement_output_file = ""; // Explicitly empty
    }

    std::string final_svg_output_file = cli_svg_output_file;
    if (final_svg_output_file.empty())
    {
        final_svg_output_file = dataset_name + "_" + effective_label_for_output + "_" + router_algo_str + "_final.svg";
    }

    std::string final_evaluation_output_file = cli_evaluation_output_file;
    if (final_evaluation_output_file.empty())
    {
        final_evaluation_output_file = dataset_name + "_" + effective_label_for_output + "_" + router_algo_str + "_eval.txt";
    }

    try
    {
        run_single_flow(
            dataset_name,
            strategy_str, // This is the placement strategy if placement is run, ignored if cli_route_from_filepath is set
            router_algo_str,
            final_placement_output_file, // Will be empty and ignored if cli_route_from_filepath is set
            final_svg_output_file,
            final_evaluation_output_file,
            cli_route_from_filepath,
            program_start_time // Pass the overall program start time for unified timing for a single run
        );
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error in processing flow: " << e.what() << std::endl;
        auto end_time_err = std::chrono::high_resolution_clock::now();
        auto duration_err = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_err - program_start_time);
        std::cout << "Total execution time (before error): " << std::fixed << std::setprecision(3) << duration_err.count() / 1000.0 << " seconds." << std::endl;
        return 1;
    }

    auto main_flow_end_time = std::chrono::high_resolution_clock::now();
    auto main_flow_duration = std::chrono::duration_cast<std::chrono::milliseconds>(main_flow_end_time - program_start_time);
    std::cout << "Total execution time for this run: " << std::fixed << std::setprecision(3) << main_flow_duration.count() / 1000.0 << " seconds." << std::endl;

    return 0;
}