#include "SteeringData.h"
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include "Sprite.h"
#include "Arrive.h"
#include <iostream>

Arrive::Arrive(Kinematic* targetInput, Sprite* characterInput)
: target(targetInput), character(characterInput) { // Use initializer list here
}

Arrive::~Arrive() {

}

void Arrive::setTargetRadius(float targetRadiusVal) {
    targetRadius = targetRadiusVal;
}
void Arrive::setSlowRadius(float slowRadiusVal) {
    slowRadius = slowRadiusVal;
}

void Arrive::execute(float timeDelta) {
    SteeringData* sd = calculateAccelerationPointer();

    if (sd != nullptr) {
        character->update(*sd, timeDelta);
        delete sd;
        sd = nullptr;
    }
}

SteeringData Arrive::calculateAcceleration() {
    return SteeringData();
}

SteeringData* Arrive::calculateAccelerationPointer() {
    SteeringData* result = new SteeringData();
    
    sf::Vector2f directionVect = target->getPosition() - character->getPosition();
    float distance = getLengthOfVector(directionVect);
    if (distance < targetRadius) {
        delete result;
        return nullptr;
    }

    float targetSpeed;
    if (distance > slowRadius) {
        targetSpeed = maxSpeed;
    } else {
        // In slow radius
        targetSpeed = maxSpeed * (distance / slowRadius);
    }

    // Desired velocity
    sf::Vector2f targetVelocity = directionVect;
    normalizeVector(targetVelocity);
    targetVelocity *= targetSpeed;

    // Assuming we have access to the current velocity of the object, which is not shown here
    sf::Vector2f resultLinear = targetVelocity - character->getVelocityVector();
    resultLinear /= timeToTarget;

    if (getLengthOfVector(resultLinear) > maxAcceleration) {
        normalizeVector(resultLinear);
        resultLinear *= maxAcceleration;
    }
    result->setLinear(resultLinear);

    result->setAngular(0); // Assuming no angular acceleration for simplicity

    return result;
}


