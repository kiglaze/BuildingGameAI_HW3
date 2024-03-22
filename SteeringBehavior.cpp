#include "SteeringBehavior.h"

float SteeringBehavior::getLengthOfVector(sf::Vector2f vectorInput) {
    // Calculate and return the length (magnitude) of the vector
    return sqrt(pow(vectorInput.x, 2) + pow(vectorInput.y, 2));
}

void SteeringBehavior::normalizeVector(sf::Vector2f& vectorInput) {
    float length = getLengthOfVector(vectorInput);
    if (length > 0) {
        vectorInput = sf::Vector2f((vectorInput.x / length), (vectorInput.y / length));
    }
}

