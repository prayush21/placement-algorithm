#ifndef PARSER_H
#define PARSER_H

#include <string>
#include "structures.h"

class Parser
{
public:
    // Static function to parse all files based on the .aux file
    static Circuit parse_ispd2011(const std::string &aux_file_path, Circuit &circuit);

private:
    // Helper functions for parsing specific file types.
    // These would typically be called by parse_ispd2011.
    // They take the file path and the Circuit object to populate.
    static void parse_nodes(const std::string &file_path, Circuit &circuit);
    static void parse_nets(const std::string &file_path, Circuit &circuit);
    static void parse_pl(const std::string &file_path, Circuit &circuit);     // Initial/Fixed placement
    static void parse_scl(const std::string &file_path, Circuit &circuit);    // Site/Core layout
    static void parse_wts(const std::string &file_path, Circuit &circuit);    // Net weights (optional)
    static void parse_shapes(const std::string &file_path, Circuit &circuit); // Fixed node shapes (optional)
    static void parse_route(const std::string &file_path, Circuit &circuit);  // Routing grid info
};

#endif // PARSER_H
