#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <unordered_set>
#include <limits>
#include <queue>

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
    std::unordered_map<int, Node*> nodeMap;

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
        if (nodeMap.find(node->identifier) == nodeMap.end()) {
            nodes.push_back(node);
            nodeMap[node->identifier] = node;
        }
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

    void loadFromCSV(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;

        while (getline(file, line)) {
            std::stringstream linestream(line);
            int fromID, toID;
            float cost;
            char delimiter;

            linestream >> fromID >> delimiter >> toID >> delimiter >> cost;

            Node* fromNode;
            Node* toNode;

            // Check if the node already exists, if not create a new node
            if (nodeMap.find(fromID) == nodeMap.end()) {
                fromNode = new Node(fromID);
                addNode(fromNode);
            } else {
                fromNode = nodeMap[fromID];
            }

            if (nodeMap.find(toID) == nodeMap.end()) {
                toNode = new Node(toID);
                addNode(toNode);
            } else {
                toNode = nodeMap[toID];
            }

            addConnection(fromNode, toNode, cost);
        }
    }



    // Dijkstra's algorithm to find the shortest path from the source node to all other nodes
    void dijkstra(int sourceId) {
        std::unordered_map<int, float> distances;
        std::unordered_map<int, int> predecessors;
        auto comp = [&distances](int lhs, int rhs) {
            return distances[lhs] > distances[rhs];
        };
        std::priority_queue<int, std::vector<int>, decltype(comp)> pq(comp);

        // Initialize distances to infinity and source distance to 0
        for (const auto& node : nodes) {
            distances[node->identifier] = std::numeric_limits<float>::infinity();
        }
        distances[sourceId] = 0.0f;

        pq.push(sourceId);

        while (!pq.empty()) {
            int currentNodeId = pq.top();
            pq.pop();
            Node* currentNode = nodeMap[currentNodeId];

            // For each neighbor of the current node
            for (const auto& conn : getConnections(*currentNode)) {
                int neighborId = conn.toNode->identifier;
                float weight = conn.cost;
                float distanceThroughU = distances[currentNodeId] + weight;
                if (distanceThroughU < distances[neighborId]) {
                    distances[neighborId] = distanceThroughU;
                    predecessors[neighborId] = currentNodeId;
                    pq.push(neighborId);
                }
            }
        }

        // Optionally: Print the shortest distances to all nodes from source
        for (const auto& dist : distances) {
            std::cout << "Shortest distance from node " << sourceId << " to node " << dist.first << " is " << dist.second << std::endl;
        }
    }

};

int main() {
    // Create a graph
    Graph graph;
    graph.loadFromCSV("subset_airport_distances_revised.csv");


/*     // Create and add nodes to the graph
    Node* node1 = new Node(1);
    Node* node2 = new Node(2);
    graph.addNode(node1);
    graph.addNode(node2);

    // Add a connection between the nodes
    graph.addConnection(node1, node2, 10.0f); */

    // Print the graph
    graph.printGraph();

    graph.generateDotFile("graph_dot_file.dot");

    graph.dijkstra(1052907);

    // The Graph destructor will delete the nodes
    return 0;
}

// Compile with: g++ -o main main.cpp
// Run with: ./main
