#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <unordered_set>
#include <limits>
#include <queue>
#include <algorithm>
#include <cmath>

// Assume Node class has an identifier
class Node {
private:
    int identifier; // Assume there is an identifier for each node for simplicity
    float x, y;     // Coordinates of the node

public:
    // Constructor to set the node's identifier and coordinates
    Node(int id, float xCoord, float yCoord) : identifier(id), x(xCoord), y(yCoord) {}

    // Getter for identifier
    int getId() const {
        return identifier;
    }

    // Getter for x
    float getX() const {
        return x;
    }

    // Setter for x
    void setX(float xCoord) {
        x = xCoord;
    }

    // Getter for y
    float getY() const {
        return y;
    }

    // Setter for y
    void setY(float yCoord) {
        y = yCoord;
    }

    // If you need to change the identifier (which is less common), add setter for identifier as well
    void setId(int id) {
        identifier = id;
    }
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

    std::unordered_map<int, Node*> getNodeMap() {
        return nodeMap;
    }

    // Function to get a list of all connections outgoing from the given node.
    std::vector<Connection> getConnections(const Node& fromNode) const {
        std::vector<Connection> result;
        for (const auto& conn : connections) {
            if (conn.fromNode->getId() == fromNode.getId()) {
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
        if (nodeMap.find(node->getId()) == nodeMap.end()) {
            nodes.push_back(node);
            nodeMap[node->getId()] = node;
        }
    }

    void printGraph() const {
        for (const auto& node : nodes) {
            std::cout << "Node " << node->getId() << " connects to: ";
            auto connections = getConnections(*node);
            for (const auto& conn : connections) {
                std::cout << conn.toNode->getId() << " (cost: " << conn.getCost() << "), ";
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
                out << node->getId() << " -> " << conn.toNode->getId();
                out << " [label=\"" << conn.getCost() << "\"];" << std::endl;
            }
        }

        out << "}" << std::endl;
    }

    void loadFromCSV(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;

        // Ingore the first line, which is the header.
        getline(file, line);

        while (getline(file, line)) {
            std::stringstream linestream(line);
            int fromID, toID;
            float cost, xCoordFrom, yCoordFrom, xCoordTo, yCoordTo;
            char delimiter;

            linestream >> fromID >> delimiter >> toID >> delimiter >> cost >> delimiter >> xCoordFrom >> delimiter >> yCoordFrom >> delimiter >> xCoordTo >> delimiter >> yCoordTo;

            Node* fromNode;
            Node* toNode;

            // Check if the node already exists, if not create a new node
            if (nodeMap.find(fromID) == nodeMap.end()) {
                fromNode = new Node(fromID, xCoordFrom, yCoordFrom);
                addNode(fromNode);
            } else {
                fromNode = nodeMap[fromID];
            }

            if (nodeMap.find(toID) == nodeMap.end()) {
                toNode = new Node(toID, xCoordTo, yCoordTo);
                addNode(toNode);
            } else {
                toNode = nodeMap[toID];
            }

            addConnection(fromNode, toNode, cost);
        }
    }


    void dijkstra(int sourceId, std::unordered_map<int, int>& predecessors) {
        std::unordered_map<int, float> distances;
        auto comp = [&distances](int lhs, int rhs) {
            return distances[lhs] > distances[rhs];
        };
        std::priority_queue<int, std::vector<int>, decltype(comp)> pq(comp);

        // Initialize distances to infinity, source distance to 0, and predecessors map
        for (const auto& node : nodes) {
            distances[node->getId()] = std::numeric_limits<float>::infinity();
            predecessors[node->getId()] = -1; // Use -1 to indicate no predecessor
        }
        distances[sourceId] = 0.0f;

        pq.push(sourceId);

        while (!pq.empty()) {
            int currentNodeId = pq.top();
            pq.pop();

            // For each neighbor of the current node
            for (const auto& conn : getConnections(*nodeMap[currentNodeId])) {
                int neighborId = conn.toNode->getId();
                float weight = conn.cost;
                float distanceThroughU = distances[currentNodeId] + weight;
                if (distanceThroughU < distances[neighborId]) {
                    distances[neighborId] = distanceThroughU;
                    predecessors[neighborId] = currentNodeId;
                    pq.push(neighborId);
                }
            }
        }

        // Distances and predecessors are now populated
        // You can use the predecessors map to reconstruct the shortest path
    }

    // Function to reconstruct the shortest path from source to target
    std::vector<int> getShortestPath(int sourceId, int targetId, const std::unordered_map<int, int>& predecessors) {
        std::vector<int> path;
        for (int at = targetId; at != -1; at = predecessors.at(at)) {
            path.push_back(at);
        }
        std::reverse(path.begin(), path.end()); // The path is constructed in reverse

        // Check if path starts with the sourceId; if not, path is unreachable
        if (path.front() != sourceId) {
            return {}; // Return an empty path to indicate no path found
        }
        
        return path;
    }

    // Convert degrees to radians.
    double degToRad(double degrees) {
        return degrees * M_PI / 180.0;
    }

    // Calculate the distance between 2 GPS co-ordinates on the Earth.
    double calculateEarthDistance(double lat1, double lon1, double lat2, double lon2) {
        // Convert latitude and longitude from degrees to radians
        lat1 = degToRad(lat1);
        lon1 = degToRad(lon1);
        lat2 = degToRad(lat2);
        lon2 = degToRad(lon2);

        // Calculate the distance
        double distance = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1)) * 3958.8;

        return distance;
    }

    float aStarHeuristic(int id1, int id2) {
        auto nodeMap = getNodeMap(); // Assuming getNodeMap() returns std::unordered_map<int, Node*>
        
        Node* node1 = nodeMap[id1];
        Node* node2 = nodeMap[id2];
        
        if (node1 == nullptr || node2 == nullptr) return std::numeric_limits<float>::infinity();
        
        // lat1, lon1, lat2, lon2
        return calculateEarthDistance(node1->getY(), node1->getX(), node2->getY(), node2->getX());
    }


    void aStar(int sourceId, int targetId, std::unordered_map<int, int>& predecessors) {
        std::unordered_map<int, float> gScore, fScore;
        // Prioritize neighboring nodes with the lowest fScore values.
        auto comp = [&fScore](int lhs, int rhs) {
            return fScore[lhs] > fScore[rhs];
        };
        std::priority_queue<int, std::vector<int>, decltype(comp)> openSet(comp);

        // Initialize gScore (cost so far) to infinity, except for the source node
        for (const auto& node : nodes) {
            gScore[node->getId()] = std::numeric_limits<float>::infinity();
            fScore[node->getId()] = std::numeric_limits<float>::infinity();
            predecessors[node->getId()] = -1; // Use -1 to indicate no predecessor
        }
        gScore[sourceId] = 0.0f;
        fScore[sourceId] = aStarHeuristic(sourceId, targetId);

        openSet.push(sourceId);

        std::unordered_set<int> openSetItems;
        while (!openSet.empty()) {
            int currentNodeId = openSet.top();
            openSet.pop();
            openSetItems.erase(currentNodeId);

            // If the goal is reached, exit the loop
            if (currentNodeId == targetId) break;

            // For each neighbor of the current node
            for (const auto& conn : getConnections(*nodeMap[currentNodeId])) {
                int neighborId = conn.toNode->getId();
                float tentative_gScore = gScore[currentNodeId] + conn.cost;

                if (tentative_gScore < gScore[neighborId]) {
                    // This path to neighbor is better than any previous one. Record it!
                    predecessors[neighborId] = currentNodeId;
                    gScore[neighborId] = tentative_gScore;
                    fScore[neighborId] = gScore[neighborId] + aStarHeuristic(neighborId, targetId);

                    if (openSetItems.find(neighborId) == openSetItems.end()) {
                        openSet.push(neighborId);
                        openSetItems.insert(neighborId);
                    }
                }
            }
        }

        // Distances and predecessors are now populated
        // Use the predecessors map to reconstruct the shortest path
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
    //graph.printGraph();

    graph.generateDotFile("graph_dot_file.dot");

    std::unordered_map<int, int> predecessors;
    int sourceNodeId = 1145501; // Example source node ID
    ////graph.dijkstra(sourceNodeId, predecessors);
    // To get and print the shortest path from source to another node, for example, node 4
    int targetNodeId = 1564301;
    graph.aStar(sourceNodeId, targetNodeId, predecessors);

    std::vector<int> path = graph.getShortestPath(sourceNodeId, targetNodeId, predecessors);
    std::cout << "Shortest path from " << sourceNodeId << " to " << targetNodeId << ": ";
    float totalCost = 0.0;
    if (!path.empty()) {
        for (size_t i = 0; i < path.size() - 1; ++i) {
            auto connections = graph.getConnections(*graph.getNodeMap()[path[i]]);
            for (const auto& conn : connections) {
                if (conn.toNode->getId() == path[i + 1]) {
                    totalCost += conn.getCost();
                    std::cout << conn.fromNode->getId() << " to " << conn.toNode->getId() << " (cost: " << conn.getCost() << "), ";
                    break;
                }
            }
        }
        std::cout << "Total cost: " << totalCost << std::endl;
    } else {
        std::cout << "No path found." << std::endl;
    }

    // The Graph destructor will delete the nodes
    return 0;
}

// Compile with: g++ -o main main.cpp
// Run with: ./main
