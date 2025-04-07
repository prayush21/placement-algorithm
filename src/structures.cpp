#include "../include/structures.h"
#include <stdexcept> // For exceptions

// --- Circuit Helper Implementations ---

Node *Circuit::add_node(const std::string &name, NodeType type, double w, double h)
{
    if (cell_map.count(name))
    {
        // Handle error: Node already exists
        // For now, let's return the existing node pointer
        // Consider logging a warning or throwing an exception based on desired behavior
        return &cell_map.at(name);
    }
    cell_map[name] = Node(name, type, w, h);
    return &cell_map.at(name);
}

Net *Circuit::add_net(const std::string &name)
{
    if (net_map.count(name))
    {
        // Handle error: Net already exists
        return &net_map.at(name);
    }
    net_map[name] = Net(name);
    return &net_map.at(name);
}

// --- Other Structure Implementations (if needed) ---
// e.g., methods for RoutingGrid initialization, etc.