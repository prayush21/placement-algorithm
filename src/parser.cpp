#include "../include/parser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>  // For exceptions
#include <filesystem> // For path manipulation (C++17)

// Helper function to get the directory path from a full file path
std::string get_dir_path(const std::string &file_path)
{
    std::filesystem::path p(file_path);
    return p.parent_path().string();
}

// Placeholder implementations for Parser methods

Circuit Parser::parse_ispd2011(const std::string &aux_file_path, Circuit &circuit)
{
    std::cout << "[Parser] Parsing ISPD 2011 benchmark from: " << aux_file_path << std::endl;
    circuit.benchmark_name = aux_file_path; // Store aux file path as identifier
    std::cout << "[Parser] Benchmark name: " << circuit.benchmark_name << std::endl;
    std::ifstream aux_file(aux_file_path);
    if (!aux_file.is_open())
    {
        throw std::runtime_error("Error: Cannot open .aux file: " + aux_file_path);
    }

    std::string line;
    std::string base_dir = get_dir_path(aux_file_path);

    // Expecting format like: RowBasedPlacement : file1.nodes file2.nets ...
    if (std::getline(aux_file, line))
    {
        std::size_t colon_pos = line.find(':');
        if (colon_pos != std::string::npos)
        {
            // Check if the line starts with the expected prefix
            std::string prefix = line.substr(0, colon_pos);
            // Trim whitespace from prefix if necessary
            prefix.erase(prefix.find_last_not_of(" \t") + 1);
            prefix.erase(0, prefix.find_first_not_of(" \t"));

            if (prefix == "RowBasedPlacement")
            {
                std::stringstream ss(line.substr(colon_pos + 1));
                std::string filename;
                while (ss >> filename)
                {
                    if (filename.empty())
                        continue;

                    std::string full_path = base_dir + "/" + filename;
                    std::cout << "[Parser] Found file: " << full_path << std::endl; // Debug output

                    // Determine file type by extension
                    if (filename.length() > 6 && filename.substr(filename.length() - 6) == ".nodes")
                    {
                        std::cout << "[Parser] Calling parse_nodes for: " << full_path << std::endl; // Debug output
                        parse_nodes(full_path, circuit);
                    }
                    else if (filename.length() > 5 && filename.substr(filename.length() - 5) == ".nets")
                    {
                        std::cout << "[Parser] Calling parse_nets for: " << full_path << std::endl; // Debug output
                        parse_nets(full_path, circuit);
                    }
                    else if (filename.length() > 3 && filename.substr(filename.length() - 3) == ".pl")
                    {
                        std::cout << "[Parser] Calling parse_pl for: " << full_path << std::endl; // Debug output
                        parse_pl(full_path, circuit);
                    }
                    else if (filename.length() > 4 && filename.substr(filename.length() - 4) == ".scl")
                    {
                        std::cout << "[Parser] Calling parse_scl for: " << full_path << std::endl; // Debug output
                        parse_scl(full_path, circuit);
                    }
                    else if (filename.length() > 4 && filename.substr(filename.length() - 4) == ".wts")
                    {
                        std::cout << "[Parser] Calling parse_wts for: " << full_path << std::endl; // Debug output
                        parse_wts(full_path, circuit);
                    }
                    else if (filename.length() > 7 && filename.substr(filename.length() - 7) == ".shapes")
                    {
                        std::cout << "[Parser] Calling parse_shapes for: " << full_path << std::endl; // Debug output
                        parse_shapes(full_path, circuit);
                    }
                    else if (filename.length() > 6 && filename.substr(filename.length() - 6) == ".route")
                    {
                        std::cout << "[Parser] Calling parse_route for: " << full_path << std::endl; // Debug output
                        parse_route(full_path, circuit);
                    }
                    else
                    {
                        std::cout << "[Parser] Warning: Unknown file type or extension for: " << filename << std::endl;
                    }
                }
            }
            else
            {
                throw std::runtime_error("Error: Unexpected format in .aux file. Expected 'RowBasedPlacement : ...', found prefix: " + prefix);
            }
        }
        else
        {
            throw std::runtime_error("Error: Colon not found in the first line of .aux file.");
        }
    }
    else
    {
        throw std::runtime_error("Error: Could not read the first line from .aux file.");
    }

    aux_file.close();
    std::cout << "[Parser] Finished parsing." << std::endl;
    return circuit;
}

// --- Private Helper Implementations (Stubs) ---
void Parser::parse_nodes(const std::string &file_path, Circuit &circuit)
{
    std::ifstream nodes_file(file_path);
    if (!nodes_file.is_open())
    {
        throw std::runtime_error("Error: Cannot open nodes file: " + file_path);
    }

    std::string line;
    while (std::getline(nodes_file, line))
    {
        if (line.empty() || line[0] == '#' ||
            line.find("Num") != std::string::npos ||
            line.find("UCLA") != std::string::npos)
        {
            continue;
        }

        // Format e.g. "cell0 width height"
        std::stringstream ss(line);
        std::string node_name;
        long long width, height;
        ss >> node_name >> width >> height;
        if (node_name.empty())
        {
            continue;
        }

        // Create node with parsed dimensions
        Node node(node_name, width, height, width * height);
        circuit.cell_map[node_name] = node;
    }
    nodes_file.close();

    std::cout << "[Parser] Parsed nodes from: " << file_path << std::endl;
}

// void Parser::parse_nodes(const std::string &file_path, Circuit &circuit)
// {
//     std::cout << "[Parser] Parsing nodes: " << file_path << " (STUB)" << std::endl;
//     // TODO: Implement .nodes file parsing logic
//     // Read lines, extract node name, width, height

//     // Handle terminal/fixed status
//     // Use circuit.add_node()
// }

void Parser::parse_nets(const std::string &file_path, Circuit &circuit)
{
    std::ifstream fin(file_path);
    if (!fin.is_open())
    {
        throw std::runtime_error("Error: Cannot open net file: " + file_path);
    }

    std::string line;
    std::string currentNetId;
    while (std::getline(fin, line))
    {
        if (line.empty() || line[0] == '#' ||
            line.find("Num") != std::string::npos ||
            line.find("UCLA") != std::string::npos)
        {
            continue;
        }

        std::stringstream ss(line);
        std::string token;
        ss >> token;

        if (token == "NetDegree")
        {
            // Format: "NetDegree : X net1234"
            std::string colon;
            int netDegree;
            ss >> colon >> netDegree >> currentNetId;

            // Create new net using Circuit's helper
            Net *net = circuit.add_net(currentNetId);
            net->partition_A_count = 0;
            net->partition_B_count = 0;
        }
        else
        {
            std::string cellId = token;
            if (circuit.cell_map.find(cellId) == circuit.cell_map.end())
            {
                // Possibly a fixed cell/pad not in cell_map
                continue;
            }

            // Add node to net's node list
            Node &cell = circuit.cell_map[cellId];
            circuit.net_map[currentNetId].nodes.push_back(&cell);

            // Update partition counts based on cell's partition
            // if (cell.partition_id == 0)
            // {
            //     circuit.net_map[currentNetId].partition_A_count++;
            // }
            // else if (cell.partition_id == -1)
            // {
            //     circuit.net_map[currentNetId].partition_B_count++;
            // }
        }
    }
    fin.close();

    std::cout << "[Parser] Read " << circuit.net_map.size()
              << " nets from " << file_path << std::endl;
}

void Parser::parse_pl(const std::string &file_path, Circuit &circuit)
{
    std::cout << "[Parser] Parsing placement: " << file_path << std::endl;
    std::ifstream pl_file(file_path);
    if (!pl_file.is_open())
    {
        throw std::runtime_error("Error: Cannot open placement file: " + file_path);
    }

    std::string line;
    while (std::getline(pl_file, line))
    {
        if (line.empty() || line[0] == '#' ||
            line.find("Num") != std::string::npos ||
            line.find("UCLA") != std::string::npos)
        {
            continue;
        }

        std::stringstream ss(line);
        std::string node_name;
        double x, y;
        std::string colon, orient, fixed_status;

        // Parse the line format: node_name x y : N /FIXED or /FIXED_NI
        ss >> node_name >> x >> y >> colon >> orient >> fixed_status;

        if (circuit.cell_map.find(node_name) == circuit.cell_map.end())
        {
            continue;
        }

        Node &cell = circuit.cell_map[node_name];
        cell.pos.x = x;
        cell.pos.y = y;
        // Set type to FIXED if either FIXED or FIXED_NI is present
        if (fixed_status == "/FIXED" || fixed_status == "/FIXED_NI")
        {
            cell.type = FIXED;
        }
        else
        {
            cell.type = MOVABLE;
        }
    }
    pl_file.close();
}

void Parser::parse_scl(const std::string &file_path, Circuit &circuit)
{
    std::cout << "[Parser] Parsing site layout: " << file_path << " (STUB)" << std::endl;
    std::ifstream scl_file(file_path);
    if (!scl_file.is_open())
    {
        throw std::runtime_error("Error: Cannot open .scl file: " + file_path);
    }
    std::string line;
    // Skip header line
    if (!std::getline(scl_file, line) || line.find("UCLA scl") == std::string::npos) // More flexible header check
    {
        throw std::runtime_error("Error: Invalid or missing .scl file header in " + file_path);
    }

    // Find NumRows
    int num_rows = 0;
    while (std::getline(scl_file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::string keyword;
        std::string colon;
        std::istringstream iss(line);
        iss >> keyword >> colon;
        if (keyword == "NumRows" && colon == ":")
        {
            if (!(iss >> num_rows))
            {
                throw std::runtime_error("Error: Failed to parse NumRows value in " + file_path);
            }
            // std::cout << "NumRows:" << num_rows << std::endl; // Remove debug print
            break; // Found NumRows, exit loop
        }
        // If we reach here without finding NumRows after non-comment lines, it's an error
        if (!keyword.empty())
        {
            throw std::runtime_error("Error: Expected 'NumRows : <count>' line not found near '" + line + "' in " + file_path);
        }
    }

    if (num_rows <= 0)
    {
        throw std::runtime_error("Error: Invalid or zero NumRows specified in " + file_path);
    }

    // std::cout << "NumRows: " << num_rows << std::endl; // Remove debug print
    circuit.core_rows.reserve(num_rows); // Pre-allocate space

    // Variables to store properties of the first and last rows for core_region calculation
    bool first_row_processed = false;
    double first_row_subrow_origin = 0.0;
    double first_row_coordinate = 0.0;
    double last_row_subrow_origin = 0.0;
    double last_row_sitewidth = 0.0;
    int last_row_num_sites = 0;
    double last_row_coordinate = 0.0;
    double last_row_height = 0.0;
    double last_row_sitespacing = 0.0;
    int current_row_index = 0;

    // Process each CoreRow
    while (std::getline(scl_file, line))
    {

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#')
            continue;

        // Trim leading/trailing whitespace (optional but good practice)
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

        if (line == "CoreRow Horizontal")
        {
            if (current_row_index >= num_rows)
            {
                std::cerr << "[Parser] Warning: Found more CoreRow blocks than specified by NumRows in " << file_path << std::endl;
                // Decide whether to break or continue parsing extra rows
                // break; // Option: Stop parsing if more rows than expected
            }

            double coordinate = 0.0;
            double height = 0.0;
            double sitewidth = 1.0;         // Default value
            double sitespacing = 1.0;       // Default value
            std::string siteorient = "N";   // Default value
            std::string sitesymmetry = "Y"; // Default value
            double subrow_origin = 0.0;
            int num_sites = 0;
            std::string end_token;

            // Read CoreRow attributes until "End"
            while (std::getline(scl_file, line) && line.find("End") == std::string::npos)
            {
                // Skip empty lines and comments within the block
                if (line.empty() || line[0] == '#')
                    continue;

                std::istringstream iss(line);
                std::string key, colon, value_str;
                iss >> key;

                if (key == "Coordinate")
                {
                    iss >> colon >> coordinate;
                }
                else if (key == "Height")
                {
                    iss >> colon >> height;
                }
                else if (key == "Sitewidth")
                {
                    iss >> colon >> sitewidth;
                }
                else if (key == "Sitespacing")
                {
                    iss >> colon >> sitespacing;
                }
                else if (key == "Siteorient")
                {
                    iss >> colon >> siteorient;
                }
                else if (key == "Sitesymmetry")
                {
                    iss >> colon >> sitesymmetry;
                }
                else if (key == "SubrowOrigin")
                {
                    // Format: SubrowOrigin : <value> NumSites : <value>
                    std::string num_sites_key, num_sites_colon;
                    iss >> colon >> subrow_origin >> num_sites_key >> num_sites_colon >> num_sites;
                    if (num_sites_key != "NumSites" || num_sites_colon != ":")
                    {
                        throw std::runtime_error("Error: Malformed SubrowOrigin/NumSites line in " + file_path + ": " + line);
                    }
                }
                else
                {
                    // Handle unknown key or malformed line if necessary
                    std::cerr << "[Parser] Warning: Unknown or malformed key '" << key << "' in CoreRow block in " << file_path << std::endl;
                }

                if (iss.fail() && !iss.eof())
                { // Check for parsing errors
                    throw std::runtime_error("Error: Failed to parse attribute line in CoreRow block: " + line);
                }
            }

            // Validate that "End" was found correctly
            if (line.find("End") == std::string::npos)
            {
                throw std::runtime_error("Error: CoreRow block not terminated with 'End' in " + file_path);
            }

            // std::cout << "Line" << current_row_index << std::endl; // Remove debug print

            // Store the first row's properties
            if (!first_row_processed)
            {
                first_row_subrow_origin = subrow_origin;
                first_row_coordinate = coordinate;
                first_row_processed = true;
            }

            // Always update the last row's properties
            last_row_subrow_origin = subrow_origin;
            last_row_sitewidth = sitewidth;
            last_row_num_sites = num_sites;
            last_row_coordinate = coordinate;
            last_row_height = height;
            last_row_sitespacing = sitespacing;

            // Store the row information in circuit data structure
            // Assuming CoreRow constructor: CoreRow(int id, double coord, double subOrigin, int nSites, double siteW, double siteSp, double h, double totalW)
            // Calculate total width if needed, or pass 0.0 if not directly stored/calculated here
            double total_width = num_sites * sitewidth + (num_sites > 0 ? (num_sites - 1) * sitespacing : 0);
            circuit.core_rows.emplace_back(current_row_index, coordinate, subrow_origin, num_sites, sitewidth, sitespacing, height, total_width);
            current_row_index++;
        }
        else if (!line.empty())
        {
            // Handle unexpected lines outside of CoreRow blocks
            std::cerr << "[Parser] Warning: Unexpected line encountered in " << file_path << ": " << line << std::endl;
        }
    }

    if (current_row_index != num_rows)
    {
        std::cerr << "[Parser] Warning: Number of CoreRows found (" << current_row_index
                  << ") does not match NumRows specified (" << num_rows << ") in " << file_path << std::endl;
    }

    // std::cout << "scl read finish: " << num_rows << std::endl; // Remove debug print

    // Calculate the core region based on the first and last rows found
    if (first_row_processed)
    {
        circuit.core_region.bottom_left.x = first_row_subrow_origin;
        circuit.core_region.bottom_left.y = first_row_coordinate;
        circuit.core_region.top_right.x = last_row_subrow_origin + last_row_sitewidth * last_row_num_sites;
        circuit.core_region.top_right.y = last_row_coordinate + last_row_height;
    }
    else
    {
        // Handle case where no rows were found, maybe set to zero or throw error
        std::cerr << "[Parser] Warning: No CoreRows found or processed in .scl file: " << file_path << ". Core region might be invalid." << std::endl;
        circuit.core_region.bottom_left.x = 0;
        circuit.core_region.bottom_left.y = 0;
        circuit.core_region.top_right.x = 0;
        circuit.core_region.top_right.y = 0;
    }

    scl_file.close(); // Good practice to close the file explicitly

    // Print the corerows:
    std::cout << "Number of corerows:" << circuit.core_rows.size() << std::endl;

    // Printing the core region
    std::cout << "Core Region Bottom Left:" << circuit.core_region.bottom_left.x << ", " << circuit.core_region.bottom_left.y << std::endl;
    std::cout << "Core Region Top Right:" << circuit.core_region.top_right.x << ", " << circuit.core_region.top_right.y << std::endl;

    // Printing the area of the code region:
    // double width = abs(circuit.core_region.bottom_left.x - circuit.core_region.top_right.x); // Use calculated width directly
    // double height = abs(circuit.core_region.bottom_left.y - circuit.core_region.top_right.y); // Use calculated height directly
    // std::cout << "[Core region Area:]" << width * height << std::endl; // Remove debug print
    std::cout << "[Parser] Parsed " << circuit.core_rows.size() << " core rows from: " << file_path << std::endl;
}

void Parser::parse_wts(const std::string &file_path, Circuit &circuit)
{
    std::cout << "[Parser] Parsing net weights: " << file_path << " (STUB)" << std::endl;
    // TODO: Implement .wts file parsing logic
    // Read net name and weight
    // Update Net.weight
}

void Parser::parse_shapes(const std::string &file_path, Circuit &circuit)
{
    std::cout << "[Parser] Parsing shapes: " << file_path << " (STUB)" << std::endl;
    // TODO: Implement .shapes file parsing logic
    // Read fixed node names and their rectangular shapes
    // Update Node.shapes for fixed nodes
}

void Parser::parse_route(const std::string &file_path, Circuit &circuit)
{
    std::cout << "[Parser] Parsing routing grid: " << file_path << " (STUB)" << std::endl;
    // TODO: Implement .route file parsing logic
    // Read grid dimensions (x, y, layers)
    // Read tile size, origin
    // Read edge capacities, min width/spacing
    // Populate Circuit.routing_grid
}