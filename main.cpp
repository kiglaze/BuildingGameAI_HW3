#include <iostream>
#include <vector>
#include <fstream>

// Assume Node class has an identifier
class Node {
public:
    int identifier; // Assume there is an identifier for each node for simplicity

    Node(int id) : identifier(id) {} // Constructor to set the node's identifier
};

// Forward declaration of Node to satisfy the Connection's need for Node
class Connection {
public:
    Node* fromNode; // The node that this connection came from.
    Node* toNode;   // The node that this connection leads to.
    float cost;     // The non-negative cost of this connection.

    // Constructor
    Connection(Node* from, Node* to, float c) : fromNode(from), toNode(to), cost(c) {}

    // Getter for the cost of the connection
    float getCost() const {
        return cost;
    }
};

class Graph {
private:
    std::vector<Node*> nodes; // A list of all nodes in the graph
    std::vector<Connection> connections; // A list of all connections

public:
    ~Graph() {
        // Make sure to delete all nodes to prevent memory leaks
        for (Node* node : nodes) {
            delete node;
        }
    }

    // Function to get a list of all connections outgoing from the given node.
    std::vector<Connection> getConnections(const Node& fromNode) const {
        std::vector<Connection> result;
        for (const auto& conn : connections) {
            if (conn.fromNode->identifier == fromNode.identifier) {
                result.push_back(conn);
            }
        }
        return result;
    }

    // Add a connection to the graph
    void addConnection(Node* from, Node* to, float cost) {
        connections.emplace_back(from, to, cost);
    }

    // Add a node to the graph
    void addNode(Node* node) {
        nodes.push_back(node);
    }

    void printGraph() const {
        for (const auto& node : nodes) {
            std::cout << "Node " << node->identifier << " connects to: ";
            auto connections = getConnections(*node);
            for (const auto& conn : connections) {
                std::cout << conn.toNode->identifier << " (cost: " << conn.getCost() << "), ";
            }
            std::cout << std::endl;
        }
    }

    void generateDotFile(const std::string& filename) {
        std::ofstream out(filename);
        out << "digraph G {" << std::endl;

        for (const auto& node : nodes) {
            auto connections = getConnections(*node);
            for (const auto& conn : connections) {
                out << node->identifier << " -> " << conn.toNode->identifier;
                out << " [label=\"" << conn.getCost() << "\"];" << std::endl;
            }
        }

        out << "}" << std::endl;
    }
};

int main() {
    // Create a graph
    Graph graph;

    // Create and add nodes to the graph
    Node* node1 = new Node(1);
    Node* node2 = new Node(2);
    graph.addNode(node1);
    graph.addNode(node2);

    // Add a connection between the nodes
    graph.addConnection(node1, node2, 10.0f);

    // Print the graph
    graph.printGraph();

    graph.generateDotFile("graph_dot_file.dot");

    // The Graph destructor will delete the nodes
    return 0;
}

// Compile with: g++ -o main main.cpp
// Run with: ./main
