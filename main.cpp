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

    void aStar(int sourceId, int targetId, std::function<float(int, int)> heuristicFun, std::unordered_map<int, int>& predecessors) {
        std::unordered_map<int, float> gScore, fScore;
        auto nodeMap = getNodeMap(); // Assume this exists and is populated elsewhere

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

        // Distances and predecessors are now populated
        // Further processing can be done outside this function
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

    void fun2(void (*fun)(int, int), int someVal1, int someVal2) {
        fun(someVal1, someVal2);
    }

    float aStarHeuristicLatLonDist(int id1, int id2) {
        auto nodeMap = getNodeMap(); // Assuming getNodeMap() returns std::unordered_map<int, Node*>
        
        Node* node1 = nodeMap[id1];
        Node* node2 = nodeMap[id2];
        
        if (node1 == nullptr || node2 == nullptr) return std::numeric_limits<float>::infinity();
        
        return calculateEarthDistance(node1->getY(), node1->getX(), node2->getY(), node2->getX());
    }

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

    void aStarLatLonDist(int sourceId, int targetId, std::unordered_map<int, int>& predecessors) {
        // Here, 'this' is valid because we are in a non-static member function
        this->aStar(sourceId, targetId,
                    [this](int id1, int id2) { return this->aStarHeuristicLatLonDist(id1, id2); }, predecessors);
    }

    void aStarEarthManhattanDist(int sourceId, int targetId, std::unordered_map<int, int>& predecessors) {
        // Here, 'this' is valid because we are in a non-static member function
        this->aStar(sourceId, targetId,
                    [this](int id1, int id2) { return this->aStarHeuristicEarthManhattanDist(id1, id2); }, predecessors);
    }

};

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

void printIt1(int id1, int id2) {
    std::cout << "I like the number " << id1 << " and the number " << id2 << std::endl;
}

// Function to add a wall
void addWall(std::vector<sf::RectangleShape>& walls, const sf::Vector2f& size, const sf::Vector2f& position, const sf::Color& color) {
    sf::RectangleShape wall(size);
    wall.setPosition(position);
    wall.setFillColor(color);
    walls.push_back(wall);
}

int convertPixelToTileNum(float pixelDimensionVal, float tileSize) {
    return (pixelDimensionVal - (tileSize / 2)) / tileSize;
}


int main()
{
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

    int sourceNodeId = 1145501; // Example source node ID
    int targetNodeId = 1564301;

    //graph.fun2(printIt1, sourceNodeId, targetNodeId);
    
    std::unordered_map<int, int> predecessorsDijkstra;
    graph.dijkstra(sourceNodeId, predecessorsDijkstra);
    printGraphShortestPath(graph, sourceNodeId, targetNodeId, predecessorsDijkstra);
    
    std::unordered_map<int, int> predecessorsAStarLatLon;
    graph.aStarLatLonDist(sourceNodeId, targetNodeId, predecessorsAStarLatLon);
    printGraphShortestPath(graph, sourceNodeId, targetNodeId, predecessorsAStarLatLon);

    std::unordered_map<int, int> predecessorsAStarManhattan;
    graph.aStarEarthManhattanDist(sourceNodeId, targetNodeId, predecessorsAStarManhattan);
    printGraphShortestPath(graph, sourceNodeId, targetNodeId, predecessorsAStarManhattan);



    // Pasting in code from HW2
        sf::Clock clock;
    float maxWindowX = 640.0;
    // Window height is 3/4 of width.
    float maxWindowY = 480.0;
    sf::RenderWindow window(sf::VideoMode(maxWindowX, maxWindowY), "SFML works!");
    
    std::vector<sf::RectangleShape> walls;
    float wallWidth = 10.0;

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


    // Getting positions of centers of tiles. Green dots mark them.
    float tileSize = 40;
    float maxTilesX = maxWindowX / tileSize;
    float maxTilesY = maxWindowY / tileSize;

    // Create a vector to hold our green dots
    std::vector<sf::CircleShape> greenDots;

    // Example positions for green dots
    std::vector<sf::Vector2f> positions = {};
/*         {100, 100}, {200, 200}, {300, 300}, {400, 400}, {500, 500}
    }; */

    for (int i = 0; i < maxTilesY; ++i) {
        for (int j = 0; j < maxTilesX; ++j) {
            bool isInLowerLeftRoom = j > convertPixelToTileNum(50, tileSize) && j <= convertPixelToTileNum(360, tileSize) && i > convertPixelToTileNum(238, tileSize) &&  i <= convertPixelToTileNum(451, tileSize);
            bool isInUpperRoom = j > convertPixelToTileNum(131, tileSize) && j <= convertPixelToTileNum(568, tileSize) && i > convertPixelToTileNum(45, tileSize) &&  i <= convertPixelToTileNum(225, tileSize);
            bool isInLowerRightRoom = j > convertPixelToTileNum(374, tileSize) && j <= convertPixelToTileNum(625, tileSize) && i > convertPixelToTileNum(238, tileSize) &&  i <= convertPixelToTileNum(451, tileSize);

            if (isInLowerLeftRoom || isInUpperRoom || isInLowerRightRoom) {
                positions.push_back(sf::Vector2f((j * tileSize) + (tileSize / 2), (i * tileSize) + (tileSize / 2)));
            }
        }
    }
/*     for (int j = 0; j < maxTilesY; ++j) {
        positions.push_back(sf::Vector2f((j * tileSize) + (tileSize / 2), (tileSize / 2)));
    } */

    // for i in 1, maxTilesX
        // push this to positions: { (tileSize/2), ((i * tileSize) - (tileSize / 2)) }

    // Populate the vector with green dots
    for (const auto& pos : positions)
    {
        sf::CircleShape dot(5); // Dot with radius 5
        dot.setFillColor(sf::Color::Green);
        dot.setPosition(pos);
        greenDots.push_back(dot);
    }



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
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);

            // Optional: Display the mouse position in the console
            std::cout << "Mouse position: " << mousePos.x << ", " << mousePos.y << std::endl;
            std::flush(std::cout); // Flush to update the position in the console in real-time

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
                // If left mouse click, arrive and align to that spot.

                if (kinemMouseClickObj != nullptr) {
                    delete kinemMouseClickObj;
                    kinemMouseClickObj = nullptr;
                }
                kinemMouseClickObj = new Kinematic(sf::Vector2f(localPosition.x, localPosition.y), 0, sf::Vector2f(0, 0), 0);
            }

            if (timeDelta > 0) {

                if (kinemMouseClickObj != nullptr) {
                    Arrive arriveBehavior(kinemMouseClickObj, spriteB);
                    arriveBehavior.execute(timeDelta);
                    // std::cout << spriteB->getVelocityVector().x << ", " << spriteB->getVelocityVector().y << std::endl;

                    Face faceBehavior(kinemMouseClickObj, spriteB);
                    faceBehavior.execute(timeDelta);

                    spriteB->dropSomeCrumbs();
                }

            }


            //myCollection.deleteMarkedSprites();
            steeringCollection.deleteMarkedSprites();
            
        }

        for(int i = 0; i < static_cast<int>(breadcrumbs.size()); i++) {
            breadcrumbs[i].draw(&window);
        }
        //myCollection.drawAll(window);
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


// Compile with: g++ -o main main.cpp
// Run with: ./main
