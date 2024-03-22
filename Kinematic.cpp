#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include "Kinematic.h"
#include "SteeringData.h"

// Constructor with member initializer list
Kinematic::Kinematic(sf::Vector2f pos, float orient, sf::Vector2f vel, float rot) {
    setPosition(pos.x, pos.y);
    setDirection(orient);
    setVelocityVector(vel.x, vel.y);
    setAngularVelocity(rot);
}
// : position(pos), direction(orient), velocity(vel), angular_velocity(rot) {}

// Get the position vector of the sprite.
sf::Vector2f Kinematic::getPosition() const {
    return position;
}

// Set the position. 0, 0 is top left of the window. down is positive y. right is positive x.
void Kinematic::setPosition(float posX, float posY) {
    position = sf::Vector2f(posX, posY);
}

// Get direction. 0 = right. 90 = down. 180 = left. 270 = up.
float Kinematic::getDirection() const {
    return direction;
}

float Kinematic::mapToRange(float rotation) {
    // Normalize the rotation to the range [-180, 180]
    rotation = fmod(rotation, 360.0f);
    if (rotation > 180.0f) {
        rotation -= 360.0f;
    } else if (rotation < -180.0f) {
        rotation += 360.0f;
    }
    return rotation;
}

// Set direction of the sprite. 0 = right. 90 = down. 180 = left. 270 = up.
void Kinematic::setDirection(float newDirection) {
    direction = mapToRange(newDirection);
}

sf::Vector2f Kinematic::getVelocityVector() const {
    return velocity;
}

void Kinematic::setVelocityVector(float valX, float valY) {
    velocity = sf::Vector2f(valX, valY);
}

float Kinematic::getAngularVelocity() const {
    return angular_velocity;
}

void Kinematic::setAngularVelocity(float w) {
    angular_velocity = w;
}

void Kinematic::update(sf::Vector2f positionVal, float orientationVal, sf::Vector2f velocityVal, float rotationVal) {
    setPosition(positionVal.x, positionVal.y);
    setDirection(orientationVal);
    setVelocityVector(velocityVal.x, velocityVal.y);
    setAngularVelocity(rotationVal);
}
void Kinematic::update(SteeringData sd, float deltaTime) {
    sf::Vector2f currentPositionVect = getPosition();
    sf::Vector2f currentVelocityVect = getVelocityVector();
    float currentOrientation = getDirection();
    float currentAngularVelocity = getAngularVelocity();

    // Position and Orientation.
    sf::Vector2f newPositionVect = currentPositionVect + (currentVelocityVect * deltaTime);
    setPosition(newPositionVect.x, newPositionVect.y);
    float newOrientation = currentOrientation + (currentAngularVelocity * deltaTime);
    setDirection(newOrientation);

    // Velocity and Rotation.
    sf::Vector2f newVelocityVect = currentVelocityVect + (sd.linear * deltaTime);
    setVelocityVector(newVelocityVect.x, newVelocityVect.y);
    float newAngularVelocity = currentAngularVelocity + (sd.angular * deltaTime);
    setAngularVelocity(newAngularVelocity);
}

