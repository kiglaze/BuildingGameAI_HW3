/* 
HW3 Assignment for CSC584: Building Game AI
Programmer: Iris Glaze
Data input csv files: full_airport_distances_revised.csv and subset_airport_distances_revised.csv
 */
#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <unordered_set>
#include <limits>
#include <queue>
#include <algorithm>
#include <cmath>
#include <typeinfo>
#include <ctime>

#include <SFML/Graphics.hpp>
#include <memory>
#include "Sprite.h"
#include "SpriteCollection.h"
#include "Arrive.h"
#include "Align.h"
#include "Face.h"
#include "SteeringData.h"
#include "Kinematic.h"
#include "LookWhereYoureGoing.h"
#include "Crumb.h"

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


// Directional connection between two Node objects.
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

// Represents a directional weighted graph which contains nodes and connections between certain nodes.
class Graph {
private:
    std::vector<Node*> nodes; // A list of all nodes in the graph
    std::vector<Connection> connections; // A list of all connections
    std::unordered_map<int, Node*> nodeMap;

public:
    // Default constructor.
    Graph() {

    }
    // Destructor
    ~Graph() {
        // Make sure to delete all nodes to prevent memory leaks
        for (Node* node : nodes) {
            delete node;
        }
    }
    // Copy constructor
    Graph(const Graph& other) {
        // Copying nodes
        for (const auto& nodePair : other.nodeMap) {
            Node* newNode = new Node(*nodePair.second); // Assuming Node has a suitable copy constructor
            this->addNode(newNode);
        }
        
        // Copying connections
        for (const auto& connection : other.connections) {
            Node* fromNode = this->nodeMap[connection.fromNode->getId()];
            Node* toNode = this->nodeMap[connection.toNode->getId()];
            this->addConnection(fromNode, toNode, connection.cost);
        }
    }

    std::unordered_map<int, Node*> getNodeMap() {
        return nodeMap;
    }

    // Add a node to the graph
    void addNode(Node* node) {
        if (nodeMap.find(node->getId()) == nodeMap.end()) {
            nodes.push_back(node);
            nodeMap[node->getId()] = node;
        }
    }

    // Add a connection to the graph
    void addConnection(Node* from, Node* to, float cost) {
        connections.emplace_back(from, to, cost);
    }

    size_t countConnections() const {
        return connections.size();
    }

    size_t countNodes() const {
        return nodes.size();
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

    // Print out the graph for debugging purposes.
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

    // Generate the .dot file to then generate the graph diagram with graphviz. 
    // Assistance obtained from ChatGPT.
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

    // Load a graph's nodes and connections from an input .csv file.
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

    // Gets the largest id of nodes that exist in the node map of this graph.
    int getLargestIdInNodeMap() {
        // Code line obtained from ChatGPT.
        auto it = std::max_element(nodeMap.begin(), nodeMap.end(),
                            [](const auto& a, const auto& b) {
                                return a.first < b.first;
                            });
        return it->first;
    }

    // Function to find a node by its x and y coordinates
    Node* findNodeByPosition(float x, float y) const {
        for (Node* node : nodes) {
            if (std::fabs(node->getX() - x) < std::numeric_limits<float>::epsilon() && 
                std::fabs(node->getY() - y) < std::numeric_limits<float>::epsilon()) {
                return node;
            }
        }
        return nullptr; // Return nullptr if no matching node is found
    }

    // Load a Graph's nodes from an input array of x, y positions.
    // The starting id offset prefents conflicting ids.
    void loadFromNodesArr(const std::vector<sf::Vector2f>& nodePositions, int startingIdOffset = 0) {
        int unusedId = 1;
        if (!nodeMap.empty()) {
            unusedId = getLargestIdInNodeMap() + 1;
        }
        unusedId = unusedId + startingIdOffset;

        for (const auto& nodePosition : nodePositions) {
            Node* additionalNode = new Node(unusedId++, nodePosition.x, nodePosition.y);
            addNode(additionalNode);
        }
    }

    // Dijkstra's Algorithm for shortest path.
    // Would need to run eithe getShortestPath() or printGraphShortestPath() afterwards to view results.
    // Some assistance from ChatGPT.
    // Results get stored in the predecessors unordered map.
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

        size_t pqSize = pq.size();
        std::cout << "# NODES IN FRINGE: " << pqSize << std::endl;

        // Distances and predecessors are now populated
        // You can use the predecessors map to reconstruct the shortest path
    }

    // A* Algorithm for shortest path. 
    // Would need to run eithe getShortestPath() or printGraphShortestPath() afterwards to view results.
    // Some assistance from ChatGPT.
    // Results get stored in the predecessors unordered map.
    void aStar(int sourceId, int targetId, std::function<float(int, int)> heuristicFun, std::unordered_map<int, int>& predecessors) {
        std::unordered_map<int, float> gScore, fScore;
        auto nodeMap = getNodeMap(); // Assume this exists and is populated elsewhere
        std::vector<Node*> fringe;

        auto comp = [&fScore](int lhs, int rhs) {
            return fScore[lhs] > fScore[rhs];
        };
        std::priority_queue<int, std::vector<int>, decltype(comp)> openSet(comp);

        for (const auto& pair : nodeMap) {
            int nodeId = pair.first;
            gScore[nodeId] = std::numeric_limits<float>::infinity();
            fScore[nodeId] = std::numeric_limits<float>::infinity();
            predecessors[nodeId] = -1; // Use -1 to indicate no predecessor
        }
        gScore[sourceId] = 0.0f;
        fScore[sourceId] = heuristicFun(sourceId, targetId);

        openSet.push(sourceId);

        std::unordered_set<int> openSetItems {sourceId};

        while (!openSet.empty()) {
            int currentNodeId = openSet.top();
            openSet.pop();
            openSetItems.erase(currentNodeId);

            fringe.push_back(nodeMap[currentNodeId]);

            if (currentNodeId == targetId) break;

            for (const auto& conn : getConnections(*nodeMap[currentNodeId])) {
                int neighborId = conn.toNode->getId();
                float tentative_gScore = gScore[currentNodeId] + conn.cost;

                if (tentative_gScore < gScore[neighborId]) {
                    predecessors[neighborId] = currentNodeId;
                    gScore[neighborId] = tentative_gScore;
                    fScore[neighborId] = gScore[neighborId] + heuristicFun(neighborId, targetId);

                    if (openSetItems.insert(neighborId).second) {
                        openSet.push(neighborId);
                    }
                }
            }
            
        }
        std::cout << "# NODES IN FRINGE: " << fringe.size() << std::endl;
        // Distances and predecessors are now populated
        // Further processing can be done outside this function
    }

    // Function to reconstruct the shortest path from source to target.
    // This gets run after either Dijkstra's or A* is run first.
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
    // Assistance from ChatGPT for the distance function.
    // lat1, lon1, lat2, lon2
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

    // Calculate the distance between two points on a 2D grid.
    double calculateEuclideanDistance(double x1, double y1, double x2, double y2) {
        return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
    }

    // Admissible heuristic for A*
    float aStarHeuristicLatLonDist(int id1, int id2) {
        auto nodeMap = getNodeMap(); // Assuming getNodeMap() returns std::unordered_map<int, Node*>
        
        Node* node1 = nodeMap[id1];
        Node* node2 = nodeMap[id2];
        
        if (node1 == nullptr || node2 == nullptr) return std::numeric_limits<float>::infinity();
        
        return calculateEarthDistance(node1->getY(), node1->getX(), node2->getY(), node2->getX());
    }

    // Inadmissible heuristic for A*
    // Calculations for the horizontal and vertical distances on a globe obtained from ChatGPT.
    float aStarHeuristicEarthManhattanDist(int id1, int id2) {
        auto nodeMap = getNodeMap(); // Assuming getNodeMap() returns std::unordered_map<int, Node*>
        
        Node* node1 = nodeMap[id1];
        Node* node2 = nodeMap[id2];
        
        if (node1 == nullptr || node2 == nullptr) return std::numeric_limits<float>::infinity();

        double lat1 = degToRad(node1->getY());
        double lon1 = degToRad(node1->getX());
        double lat2 = degToRad(node2->getY());
        double lon2 = degToRad(node2->getX());

        double avgLat = (lat1 + lat2)/2;
        double deltaLon = fabsf(lon1 - lon2);
        double deltaLat = fabsf(lat1 - lat2);

        double earthRadiusMiles = 3958.8;

        double diffVaryingLats = earthRadiusMiles * deltaLat;
        double diffVaryingLons = earthRadiusMiles * deltaLon * cos(avgLat);

        return diffVaryingLats + diffVaryingLons;        
    }

    // Heuristic is for the SFML portion, not the flight distance data.
    float aStarHeuristicEuclideanDist(int id1, int id2) {
        auto nodeMap = getNodeMap(); // Assuming getNodeMap() returns std::unordered_map<int, Node*>
        
        Node* node1 = nodeMap[id1];
        Node* node2 = nodeMap[id2];
        
        if (node1 == nullptr || node2 == nullptr) return std::numeric_limits<float>::infinity();
        
        return calculateEuclideanDistance(node1->getX(), node1->getY(), node2->getX(), node2->getY());
    }

    // Run A* with the distance between two GPS coordinates on the globe as the heuristic.
    void aStarLatLonDist(int sourceId, int targetId, std::unordered_map<int, int>& predecessors) {
        // Here, 'this' is valid because we are in a non-static member function
        this->aStar(sourceId, targetId,
                    [this](int id1, int id2) { return this->aStarHeuristicLatLonDist(id1, id2); }, predecessors);
    }

    // Run A* with the horizontal distance plus the vertical distance between two GPS coordinates on a globe.
    void aStarEarthManhattanDist(int sourceId, int targetId, std::unordered_map<int, int>& predecessors) {
        // Here, 'this' is valid because we are in a non-static member function
        this->aStar(sourceId, targetId,
                    [this](int id1, int id2) { return this->aStarHeuristicEarthManhattanDist(id1, id2); }, predecessors);
    }

    // Run A* with the standard 2D Euclidean Distance as the heuristic.
    void aStarEuclideanDist(int sourceId, int targetId, std::unordered_map<int, int>& predecessors) {
        // Here, 'this' is valid because we are in a non-static member function
        this->aStar(sourceId, targetId,
                    [this](int id1, int id2) { return this->aStarHeuristicEuclideanDist(id1, id2); }, predecessors);
    }

    // Functions for making it work with SFML.
    // Function to calculate Euclidean distance between two nodes
    float calculateDistanceBetweenNodes(const Node* a, const Node* b) {
        float dx = a->getX() - b->getX();
        float dy = a->getY() - b->getY();
        return sqrt(dx * dx + dy * dy);
    }
    // Make a graph with only nodes and no edges be fully connected, given the existing nodes.
    void connectAllNodes() {
        // Loop through all pairs of nodes
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (size_t j = i + 1; j < nodes.size(); ++j) {
                // Calculate the distance between nodes[i] and nodes[j]
                float distance = calculateDistanceBetweenNodes(nodes[i], nodes[j]);

                // Add connection between nodes[i] and nodes[j] in both directions with the calculated distance as the cost
                addConnection(nodes[i], nodes[j], distance);
                addConnection(nodes[j], nodes[i], distance);
            }
        }
    }

    // Merge another graph into this one. 
    // Assumes that the two node graphs are completely exclusive from each other.
    // The point is to add connections between the two after. Used for adding division scheme regions.
    void mergeGraph(const Graph& other) {
        int unusedId = 1;
        if (!nodeMap.empty()) {
            unusedId = getLargestIdInNodeMap() + 1;
        }
        // Iterate over all nodes in the other graph
        for (const Node* otherNode : other.nodes) {
            // Check if the node already exists in this graph
            Node* copiedNode = new Node(unusedId++, otherNode->getX(), otherNode->getY()); // Assumes Node has a copy constructor
            this->addNode(copiedNode);
            // If the node already exists, no need to add it again
        }

        // Iterate over all connections in the other graph
        for (const Connection& otherConnection : other.connections) {
            Node* fromNode = this->nodeMap[otherConnection.fromNode->getId()];
            Node* toNode = this->nodeMap[otherConnection.toNode->getId()];

            // Add connection if it doesn't already exist
            // This naive implementation adds the connection directly.
            // You may want to check for existing connections between these nodes to avoid duplicates.
            this->addConnection(fromNode, toNode, otherConnection.cost);
        }
    }

    // Allows the combining of two division schemes.
    // Must specify the exact x, y coordinates.
    void addConnection2DByCoordinates(float x1, float y1, float x2, float y2) {
        Node* n1 = findNodeByPosition(x1, y1);
        Node* n2 = findNodeByPosition(x2, y2);
        if(n1 != nullptr && n2 != nullptr) {
            float weightDistance = calculateEuclideanDistance(x1, y1, x2, y2);
            addConnection(n1, n2, weightDistance);
            addConnection(n2, n1, weightDistance);
        }
    }

    // Creates a subgraph out of the specifies locations, makes them all interconnected.
    // Then, appends the edges and nodes to the original graph.
    void performGraphAppend(std::vector<sf::Vector2f> &positionsTopRoomA) {
        Graph gameSubGraphTopRoomA;
        int largestIdGameGraph = getLargestIdInNodeMap();
        gameSubGraphTopRoomA.loadFromNodesArr(positionsTopRoomA, largestIdGameGraph);
        gameSubGraphTopRoomA.connectAllNodes();
        mergeGraph(gameSubGraphTopRoomA);
    }

};

// Given a graph that has had either Dijkstra or A* run on it, print out the shortest path and relevant costs.
void printGraphShortestPath(Graph &graph, int sourceNodeId, int targetNodeId, std::unordered_map<int, int> &predecessors)
{
    std::vector<int> path = graph.getShortestPath(sourceNodeId, targetNodeId, predecessors);
    std::cout << "Shortest path from " << sourceNodeId << " to " << targetNodeId << ": ";
    float totalCost = 0.0;
    if (!path.empty())
    {
        for (size_t i = 0; i < path.size() - 1; ++i)
        {
            auto connections = graph.getConnections(*graph.getNodeMap()[path[i]]);
            for (const auto &conn : connections)
            {
                if (conn.toNode->getId() == path[i + 1])
                {
                    totalCost += conn.getCost();
                    std::cout << conn.fromNode->getId() << " to " << conn.toNode->getId() << " (cost: " << conn.getCost() << "), ";
                    break;
                }
            }
        }
        std::cout << "Total cost: " << totalCost << std::endl;
    }
    else
    {
        std::cout << "No path found." << std::endl;
    }
}

// Function to add a wall to the SFML game.
void addWall(std::vector<sf::RectangleShape>& walls, const sf::Vector2f& size, const sf::Vector2f& position, const sf::Color& color) {
    sf::RectangleShape wall(size);
    wall.setPosition(position);
    wall.setFillColor(color);
    walls.push_back(wall);
}

// Convert the pixel in either the x or y dimension to the corresponding tile index.
int convertPixelToTileNum(float pixelDimensionVal, float tileSize) {
    return floor(pixelDimensionVal / tileSize);
}

// Convert the tile index (in either the x or y dimension) to the corresponding tile dot dimension value.
float convertTileNumToPixel(int tileNum, float tileSize) {
    return (float(tileNum) * tileSize) + (tileSize / 2);
}

// Convert the pixel (in either x or y direction) to the nearest tile dot pixel.
float convertPixelDimToTileDotDim(float pixelDim, float tileSize) {
    return convertTileNumToPixel(convertPixelToTileNum(pixelDim, tileSize), tileSize);
}

// Given the start and end times, print the time elapsed. Used for execution time comparisons.
void printTimeElapsed(clock_t timeStart, clock_t timeEnd) {
    double elapsed = double(timeEnd - timeStart) / CLOCKS_PER_SEC;
    std::cout << "Time taken by function: " 
         << elapsed << " seconds" << std::endl;
}

void performGraphAppend(int &largestIdGameGraph, std::vector<sf::Vector2f> &positionsTopRoomA, Graph &gameGraph);

// Runs multple algorithms in one go, for the sake of comparison.
void runAllShortestPathAlgs(int sourceNodeId, int targetNodeId, Graph& graph) {
    std::cout << "Running for nodes: SOURCE: " << sourceNodeId << "; TARGET: " << targetNodeId << std::endl;
    
    std::cout << "Dijkstra's" << std::endl;
    std::unordered_map<int, int> predecessorsDijkstra;
    clock_t startD = clock();
    graph.dijkstra(sourceNodeId, predecessorsDijkstra);
    clock_t stopD = clock();
    printTimeElapsed(startD, stopD);
    printGraphShortestPath(graph, sourceNodeId, targetNodeId, predecessorsDijkstra);
    std::cout << "" << std::endl;

    std::cout << "A* Admissible" << std::endl;
    std::unordered_map<int, int> predecessorsAStarLatLon;
    clock_t startA1 = clock();
    graph.aStarLatLonDist(sourceNodeId, targetNodeId, predecessorsAStarLatLon);
    clock_t stopA1 = clock();
    printTimeElapsed(startA1, stopA1);
    printGraphShortestPath(graph, sourceNodeId, targetNodeId, predecessorsAStarLatLon);
    std::cout << "" << std::endl;

    std::cout << "A* In-admissible" << std::endl;
    std::unordered_map<int, int> predecessorsAStarManhattan;
    clock_t startA2 = clock();
    graph.aStarEarthManhattanDist(sourceNodeId, targetNodeId, predecessorsAStarManhattan);
    clock_t stopA2 = clock();
    printTimeElapsed(startA2, stopA2);
    printGraphShortestPath(graph, sourceNodeId, targetNodeId, predecessorsAStarManhattan);
    std::cout << "" << std::endl;
}

int main()
{
    // Create the smaller graph.
    Graph graph;
    graph.loadFromCSV("subset_airport_distances_revised.csv");
    graph.generateDotFile("graph_dot_file.dot");

    // Create the big graph.
    Graph bigGraph;
    bigGraph.loadFromCSV("full_airport_distances_revised.csv");

    // Assess the different algorithms for the different combinations of source and target nodes.
    // Measures determined cost, execution time, and fringe node count.
    int sourceNodeId = 1145501; // Example source node ID
    int targetNodeId = 1564301;

    std::cout << "SMALL GRAPH: " << std::endl;
    runAllShortestPathAlgs(sourceNodeId, targetNodeId, graph);
    std::cout << "" << std::endl;
    std::cout << "LARGE GRAPH: " << std::endl;
    runAllShortestPathAlgs(sourceNodeId, targetNodeId, bigGraph);

    int sourceNodeId2 = 1252306;
    int targetNodeId2 = 1311102;
    std::cout << "SMALL GRAPH: " << std::endl;
    runAllShortestPathAlgs(sourceNodeId2, targetNodeId2, graph);
    std::cout << "" << std::endl;
    std::cout << "LARGE GRAPH: " << std::endl;
    runAllShortestPathAlgs(sourceNodeId2, targetNodeId2, bigGraph);

    int sourceNodeId3 = 1599403;
    int targetNodeId3 = 1067302;
    std::cout << "SMALL GRAPH: " << std::endl;
    runAllShortestPathAlgs(sourceNodeId3, targetNodeId3, graph);
    std::cout << "" << std::endl;
    std::cout << "LARGE GRAPH: " << std::endl;
    runAllShortestPathAlgs(sourceNodeId3, targetNodeId3, bigGraph);


    // START OF SFML CODE PORTION.
    sf::Clock clock;
    float maxWindowX = 640.0;
    // Window height is 3/4 of width.
    float maxWindowY = 480.0;
    sf::RenderWindow window(sf::VideoMode(maxWindowX, maxWindowY), "SFML works!");
    
    float tileSize = 40;

    std::vector<sf::RectangleShape> walls;
    float wallWidth = 10.0;

    // Add a bunch of walls to create the 3 different rooms.
    // bottom-left room: x range: 50-360. y range: 238-451
    addWall(walls, sf::Vector2f(maxWindowX * 0.55, wallWidth), sf::Vector2f(40, maxWindowY - (wallWidth + 20)), sf::Color::Black);
    addWall(walls, sf::Vector2f(wallWidth, maxWindowX * 0.35), sf::Vector2f(40, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Black);
    //addWall(walls, sf::Vector2f(maxWindowX * 0.55, wallWidth), sf::Vector2f(10, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Blue);
    addWall(walls, sf::Vector2f(maxWindowX * 0.30, wallWidth), sf::Vector2f(40, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Green);
    addWall(walls, sf::Vector2f(maxWindowX * 0.20, wallWidth), sf::Vector2f(320, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Green);
    //addWall(walls, sf::Vector2f(wallWidth, (maxWindowX * 0.35) + wallWidth), sf::Vector2f(10 + (maxWindowX * 0.55), maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Blue);
    addWall(walls, sf::Vector2f(wallWidth, (maxWindowX * 0.12) + wallWidth), sf::Vector2f(40 + (maxWindowX * 0.55), maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Green);
    addWall(walls, sf::Vector2f(wallWidth, (maxWindowX * 0.10) + wallWidth), sf::Vector2f(40 + (maxWindowX * 0.55), 385), sf::Color::Green);

    // top room: x range: 131(left)-568(right). y range: 45(top)-225(bottom)
    //addWall(walls, sf::Vector2f(maxWindowX * 0.70, wallWidth), sf::Vector2f(120, maxWindowY - (wallWidth + 20) - 225), sf::Color::Blue);
    addWall(walls, sf::Vector2f(wallWidth, maxWindowX * 0.30), sf::Vector2f(120, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.30) - 225), sf::Color::Black);
    addWall(walls, sf::Vector2f(maxWindowX * 0.70, wallWidth), sf::Vector2f(120, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.30) - 225), sf::Color::Black);
    addWall(walls, sf::Vector2f(wallWidth, (maxWindowX * 0.30) + wallWidth), sf::Vector2f(120 + (maxWindowX * 0.70), maxWindowY - (wallWidth + 20) - (maxWindowX * 0.30) - 225), sf::Color::Black);

    // bottom-right room: x range: 374-625. y range: 238-450
    addWall(walls, sf::Vector2f(maxWindowX * 0.35, wallWidth), sf::Vector2f(400, maxWindowY - (wallWidth + 20)), sf::Color::Black);
    //addWall(walls, sf::Vector2f(wallWidth, maxWindowX * 0.35), sf::Vector2f(370, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Blue);
    //addWall(walls, sf::Vector2f(maxWindowX * 0.40, wallWidth), sf::Vector2f(370, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Blue);
    addWall(walls, sf::Vector2f(maxWindowX * 0.10, wallWidth), sf::Vector2f(400, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Green);
    addWall(walls, sf::Vector2f(maxWindowX * 0.10, wallWidth), sf::Vector2f(560, maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Green);
    addWall(walls, sf::Vector2f(wallWidth, (maxWindowX * 0.35) + wallWidth), sf::Vector2f(400 + (maxWindowX * 0.35), maxWindowY - (wallWidth + 20) - (maxWindowX * 0.35)), sf::Color::Black);

    // Obstacles
    addWall(walls, sf::Vector2f(tileSize * 3.25, tileSize * 3), sf::Vector2f(50, 340), sf::Color::Black);
    addWall(walls, sf::Vector2f(tileSize * 2, tileSize), sf::Vector2f(300, 100), sf::Color::Black);
    addWall(walls, sf::Vector2f(tileSize * 2.75, tileSize), sf::Vector2f(460, 100), sf::Color::Black);


    // Getting positions of centers of tiles. Green dots mark them.
    // Tile size is defined earlier as 40.
    float maxTilesX = maxWindowX / tileSize;
    float maxTilesY = maxWindowY / tileSize;

    // Create a vector to hold our green dots
    std::vector<sf::CircleShape> greenDots;

    // Example positions for green dots
    std::vector<sf::Vector2f> positions = {};
/*         {100, 100}, {200, 200}, {300, 300}, {400, 400}, {500, 500}
    }; */
    std::vector<sf::Vector2f> positionsBottomLeft = {};
    std::vector<sf::Vector2f> positionsBottomLeftA = {};
    std::vector<sf::Vector2f> positionsBottomLeftB = {};

    std::vector<sf::Vector2f> positionsBottomRight = {};

    std::vector<sf::Vector2f> positionsTopRoom = {};
    std::vector<sf::Vector2f> positionsTopRoomA = {};
    std::vector<sf::Vector2f> positionsTopRoomB = {};
    std::vector<sf::Vector2f> positionsTopRoomC = {};
    std::vector<sf::Vector2f> positionsTopRoomD = {};

    // Here, j is x and i is y.
    for (int i = 0; i < maxTilesY; ++i) {
        for (int j = 0; j < maxTilesX; ++j) {
            bool isInLowerLeftRoom = j >= convertPixelToTileNum(50, tileSize) && j <= convertPixelToTileNum(360, tileSize) && i > convertPixelToTileNum(238, tileSize) &&  i < convertPixelToTileNum(451, tileSize);
            bool isInLowerRightRoom = j > convertPixelToTileNum(374, tileSize) && j <= convertPixelToTileNum(625, tileSize) && i > convertPixelToTileNum(238, tileSize) &&  i < convertPixelToTileNum(451, tileSize);
            bool isInUpperRoom = j >= convertPixelToTileNum(131, tileSize) && j < convertPixelToTileNum(568, tileSize) && i >= convertPixelToTileNum(45, tileSize) &&  i <= convertPixelToTileNum(225, tileSize);
            bool isInTopRoomA = isInUpperRoom && j < convertPixelToTileNum(300, tileSize);
            bool isInTopRoomB = isInUpperRoom && j >= convertPixelToTileNum(300, tileSize) && i > convertPixelToTileNum(140, tileSize);
            bool isInTopRoomC = isInUpperRoom && j >= convertPixelToTileNum(300, tileSize) && i < convertPixelToTileNum(100, tileSize);
            bool isInTopRoomD = isInUpperRoom && j == convertPixelToTileNum(420, tileSize) && i >= convertPixelToTileNum(100, tileSize) && i <= convertPixelToTileNum(140, tileSize);

            bool isInBottomLeftA = isInLowerLeftRoom && j > convertPixelToTileNum(180, tileSize);
            bool isInBottomLeftB = isInLowerLeftRoom && j <= convertPixelToTileNum(180, tileSize) && i < convertPixelToTileNum(340, tileSize);

            //sf::Vector2f dotPosVect = sf::Vector2f((j * tileSize) + (tileSize / 2), (i * tileSize) + (tileSize / 2));
            
            sf::Vector2f dotPosVect = sf::Vector2f(convertTileNumToPixel(j, tileSize), convertTileNumToPixel(i, tileSize));
            if (isInLowerLeftRoom) {
                positions.push_back(dotPosVect);
                positionsBottomLeft.push_back(dotPosVect);
            }
            if (isInLowerRightRoom) {
                positions.push_back(dotPosVect);
                positionsBottomRight.push_back(dotPosVect);
            }

            if(isInTopRoomA) {
                positions.push_back(dotPosVect);
                positionsTopRoomA.push_back(dotPosVect);
            }
            
            if(isInTopRoomB) {
                positions.push_back(dotPosVect);
                positionsTopRoomB.push_back(dotPosVect);
            }
            
            if(isInTopRoomC) {
                positions.push_back(dotPosVect);
                positionsTopRoomC.push_back(dotPosVect);
            }

            if(isInTopRoomD) {
                positions.push_back(dotPosVect);
                positionsTopRoomD.push_back(dotPosVect);
            }

            if(isInBottomLeftA) {
                positions.push_back(dotPosVect);
                positionsBottomLeftA.push_back(dotPosVect);
            }

            if(isInBottomLeftB) {
                positions.push_back(dotPosVect);
                positionsBottomLeftB.push_back(dotPosVect);
            }
        }
    }

    // Build the graph for game path finding purposes here.
    Graph gameGraph;

    // Combinining the multiple division scheme regions for the bottom left room.
    gameGraph.loadFromNodesArr(positionsBottomLeftA);
    gameGraph.connectAllNodes();

    gameGraph.performGraphAppend(positionsBottomLeftB);

    // Combinining the single regions for the bottom right room with the existing game's node graph.
    gameGraph.performGraphAppend(positionsBottomRight);

    gameGraph.addConnection2DByCoordinates(380, 340, 420, 340);

    // Combinining the multiple division scheme regions for the top room with the existing game's node graph.
    gameGraph.performGraphAppend(positionsTopRoomA);
    gameGraph.performGraphAppend(positionsTopRoomB);
    gameGraph.performGraphAppend(positionsTopRoomC);
    gameGraph.performGraphAppend(positionsTopRoomD);

    // Adding connections between certain nodes in the various division scheme graph sub-regions.
    gameGraph.addConnection2DByCoordinates(260, 180, 260, 260);
    gameGraph.addConnection2DByCoordinates(500, 180, 500, 260);

    gameGraph.addConnection2DByCoordinates(260, 180, 300, 180);
    gameGraph.addConnection2DByCoordinates(260, 260, 300, 180);

    gameGraph.addConnection2DByCoordinates(260, 60, 300, 60);
    gameGraph.addConnection2DByCoordinates(260, 100, 300, 60);

    gameGraph.addConnection2DByCoordinates(420, 180, 420, 140);
    gameGraph.addConnection2DByCoordinates(420, 100, 420, 60);

    gameGraph.addConnection2DByCoordinates(180, 260, 220, 260);
    gameGraph.addConnection2DByCoordinates(180, 300, 220, 300);
    //gameGraph.addConnection2DByCoordinates(180, 340, 220, 340);
    gameGraph.addConnection2DByCoordinates(180, 300, 220, 340);
    //gameGraph.addConnection2DByCoordinates(180, 340, 220, 380);

    // Populate the vector with green dots for visual aid purposes and debugging help.
    for (const auto& pos : positions)
    {
        sf::CircleShape dot(1); // Dot with radius 5
        dot.setFillColor(sf::Color::Green);
        dot.setPosition(pos);
        greenDots.push_back(dot);
    }

    // Part of the adding breadcrumbs process.
    std::vector<Crumb> breadcrumbs = std::vector<Crumb>();
    for(int i = 0; i < 100; i++)
    {
        Crumb c(i);
        breadcrumbs.push_back(c);
    }

    // Variables to store the mouse positions
    sf::Vector2i previousMousePositionInt = sf::Mouse::getPosition(window);
    sf::Vector2f previousMousePosition = sf::Vector2f(static_cast<float>(previousMousePositionInt.x), static_cast<float>(previousMousePositionInt.y));
    sf::Vector2f currentMousePosition;

    //SpriteCollection myCollection;
    SpriteCollection steeringCollection(&breadcrumbs);

    ////float numPixels = 1.0;
    // Frame counter
    int frameCounter = 0;
    std::string textureFilePath = "./sprite_high_res.png";

    Sprite* spriteB = new Sprite(textureFilePath, 275.f, 325.f, 0, sf::Vector2f(0, 0), 0, 0, &breadcrumbs);
    
    std::vector<Kinematic*> kinematics;
    std::vector<Sprite*> sprites;

    sprites.push_back(spriteB);
    for (Sprite* sprite : sprites) {
        kinematics.push_back(sprite);
    }

    steeringCollection.addSprite(spriteB);

    Kinematic* kinemMouseClickObj = nullptr;
    Node* nodeNearClick = nullptr;
    int* currentArriveGoalNodeId = nullptr;
    Node *currentArriveGoalNode;
    std::vector<int> gameNodeIdPath = {};
    int* gameNodeIdPathEnd = nullptr;
    std::unordered_map<int, Node*> gameGraphNodeMap = gameGraph.getNodeMap();
    Kinematic* kinemArriveGoalObj = nullptr;
    
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Increment frame counter
        frameCounter++;

        window.clear(sf::Color(200, 0, 0, 255));
        

        int frameCountMark = 10;

        // Check if 100 frames have passed
        if (frameCounter >= frameCountMark)
        {
            // Get the current mouse position in the window
            //sf::Vector2i mousePos = sf::Mouse::getPosition(window);

            // Optional: Display the mouse position in the console
            //std::cout << "\rMouse position: " << mousePos.x << ", " << mousePos.y << "      "<< std::flush;

            sf::Time elapsed = clock.restart();
            float timeDelta = elapsed.asMilliseconds();
            if (timeDelta > 0) { // Check to prevent division by zero
                // Update current mouse position
                sf::Vector2i currentMousePosInt = sf::Mouse::getPosition(window);
                currentMousePosition = sf::Vector2f(static_cast<float>(currentMousePosInt.x), static_cast<float>(currentMousePosInt.y));

                // Update previous mouse position
                previousMousePosition = currentMousePosition;

            }

            // Reset frame counter
            frameCounter = 0;

            // Activated by pressing "V" to have sprite velocity match the mouse.
            if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                sf::Vector2i localPosition = sf::Mouse::getPosition(window);
                float localPositionX = float(localPosition.x);
                float localPositionY = float(localPosition.y);
                
                // If left mouse click, arrive and align to that spot.

                if (kinemMouseClickObj != nullptr) {
                    delete kinemMouseClickObj;
                    kinemMouseClickObj = nullptr;
                }

                // Gets location of grid tile dot closest to the mouse click.
                kinemMouseClickObj = new Kinematic(sf::Vector2f(convertPixelDimToTileDotDim(localPositionX, tileSize), convertPixelDimToTileDotDim(localPositionY, tileSize)), 0, sf::Vector2f(0, 0), 0);
                sf::Vector2f kinemMouseClickObjPos = kinemMouseClickObj->getPosition();

                nodeNearClick = gameGraph.findNodeByPosition(kinemMouseClickObjPos.x, kinemMouseClickObjPos.y);

                if (nodeNearClick != nullptr) {
                    if(currentArriveGoalNodeId == nullptr) {
                        sf::Vector2f spriteBPos = spriteB->getPosition();
                        sf::Vector2f startingTileDot = sf::Vector2f(convertPixelDimToTileDotDim(spriteBPos.x, tileSize), convertPixelDimToTileDotDim(spriteBPos.y, tileSize));
                        Node* nodeNearStart = gameGraph.findNodeByPosition(startingTileDot.x, startingTileDot.y);
                        // Finding the shortest path of sprite to the click location using A* with the Euclidean distance heuristic.
                        std::unordered_map<int, int> predecessorsAStarEuclidean;
                        int startNodeId = nodeNearStart->getId();
                        int goalNodeId = nodeNearClick->getId();
                        gameGraph.aStarEuclideanDist(startNodeId, goalNodeId, predecessorsAStarEuclidean);
                        //printGraphShortestPath(gameGraph, startNodeId, goalNodeId, predecessorsAStarEuclidean);
                        gameNodeIdPath = gameGraph.getShortestPath(startNodeId, goalNodeId, predecessorsAStarEuclidean);

                        // currentArriveGoalNodeId determines the next node the sprite should focus on arriving at.
                        // It changes when a sprite reaches a node within the path that isn't the last node required to reach the mouse click.
                        currentArriveGoalNodeId = gameNodeIdPath.data();
                        gameNodeIdPathEnd = gameNodeIdPath.data() + gameNodeIdPath.size() - 1;

                        currentArriveGoalNode = gameGraphNodeMap[*currentArriveGoalNodeId];
                        if (kinemArriveGoalObj != nullptr) {
                            delete kinemArriveGoalObj;
                            kinemArriveGoalObj = nullptr; // To prevent dangling pointers
                        }
                        kinemArriveGoalObj = new Kinematic(sf::Vector2f(currentArriveGoalNode->getX(), currentArriveGoalNode->getY()), 0, sf::Vector2f(0, 0), 0);

                    }
                    
                } else {
                    if (kinemMouseClickObj != nullptr) {
                        delete kinemMouseClickObj;
                        kinemMouseClickObj = nullptr; // To prevent dangling pointers
                    }
                }

            }

            if (timeDelta > 0) {
                
                if (kinemArriveGoalObj != nullptr) {
                    Arrive arriveBehavior(kinemArriveGoalObj, spriteB);
                    arriveBehavior.execute(timeDelta);

                    Face faceBehavior(kinemArriveGoalObj, spriteB);
                    faceBehavior.execute(timeDelta);

                    spriteB->dropSomeCrumbs();
                
                    // Need to only do this once currentArriveGoalNode is reached by spriteB.
                    if(spriteB->hasArrivedAtKinemObj(kinemArriveGoalObj, 30)) {
                        // If this is not the last node in the path, set the next current goal node.
                        if(currentArriveGoalNodeId != gameNodeIdPathEnd) {
                            currentArriveGoalNodeId++;
                            currentArriveGoalNode = gameGraphNodeMap[*currentArriveGoalNodeId];
                            if (kinemArriveGoalObj != nullptr) {
                                delete kinemArriveGoalObj;
                                kinemArriveGoalObj = nullptr;
                            }
                            kinemArriveGoalObj = new Kinematic(sf::Vector2f(currentArriveGoalNode->getX(), currentArriveGoalNode->getY()), 0, sf::Vector2f(0, 0), 0);
                        // If this is the final destination node in a path, then do some memory cleanup.
                        } else {
                            if (currentArriveGoalNodeId != nullptr) {
                                currentArriveGoalNodeId = nullptr;
                            }
                            if (kinemArriveGoalObj != nullptr) {
                                delete kinemArriveGoalObj;
                                kinemArriveGoalObj = nullptr;
                            }
                            if (currentArriveGoalNode != nullptr) {
                                currentArriveGoalNode = nullptr; // To prevent dangling pointers
                            }
                            gameNodeIdPath = {};
                            
                        }
                    }
                    
                }

            }

            steeringCollection.deleteMarkedSprites();
            
        }

        for(int i = 0; i < static_cast<int>(breadcrumbs.size()); i++) {
            breadcrumbs[i].draw(&window);
        }
        steeringCollection.drawAll(window);
        // Draw all the walls
        for (const auto& wall : walls) {
            window.draw(wall);
        }

        // Draw all green dots to mark the tiles.
        for (const auto& dot : greenDots)
        {
            window.draw(dot);
        }

        window.display();
    }

    if (kinemMouseClickObj != nullptr) {
        delete kinemMouseClickObj;
        kinemMouseClickObj = nullptr; // To prevent dangling pointers
    }


    // The Graph destructor will delete the nodes
    return 0;
}

