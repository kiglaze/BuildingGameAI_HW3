// Include the header file
#include "SteeringData.h"

// Implement the getter for the linear vector
sf::Vector2f SteeringData::getLinear() {
    return linear;
}

// Implement the getter for the angular velocity/direction
float SteeringData::getAngular() {
    return angular;
}

// Implement the setter for the linear vector
void SteeringData::setLinear(sf::Vector2f lin) {
    linear = lin;
}

// Implement the setter for the angular velocity/direction
void SteeringData::setAngular(float ang) {
    angular = ang;
}
