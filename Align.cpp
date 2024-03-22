#include "SteeringData.h"
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include "Sprite.h"
#include "Align.h"
#include <iostream>

Align::Align(Kinematic* targetInput, Sprite* characterInput) {
    target = targetInput;
    character = characterInput;
}

Align::~Align() {

}

void Align::setTargetRadius(float targetRad) {
    targetRadius = targetRad;
}

void Align::setMaxAngularAcceleration(float maxAngAcc) {
    maxAngularAcceleration = maxAngAcc;
}
void Align::setMaxRotation(float maxRot) {
    maxRotation = maxRot;
}
void Align::setSlowRadius(float slowRadiusVal) {
    slowRadius = slowRadiusVal;
}
void Align::setTimeToTarget(float timeToTargetVal) {
    timeToTarget = timeToTargetVal;
}

float Align::mapToRange(float rotation) {
    // Normalize the rotation to the range [-180, 180]
    rotation = fmod(rotation, 360.0f);
    if (rotation > 180.0f) {
        rotation -= 360.0f;
    } else if (rotation < -180.0f) {
        rotation += 360.0f;
    }
    return rotation;
}

void Align::execute(float timeDelta) {
    SteeringData* sd = calculateAccelerationPointer();

    if (sd != nullptr) {
        character->update(*sd, timeDelta);
        delete sd;
        sd = nullptr;
    } else {
        character->setAngularVelocity(0.0f);
    }
}

SteeringData Align::calculateAcceleration() {
    return SteeringData();
}

SteeringData* Align::calculateAccelerationPointer() {
    SteeringData* result = new SteeringData();

    float rotation = target->getDirection() - character->getDirection();
    rotation = mapToRange(rotation);
    float rotationSize = fabs(rotation);
    
    if (rotationSize < targetRadius) {
        //std::cout << "A" << std::endl;
        delete result;
        return nullptr;
    } else {
        float targetRotation;
        if (rotationSize > slowRadius) {
            //std::cout << "B" << std::endl;
            targetRotation = maxRotation;
        } else {
            //std::cout << "C" << std::endl;
            targetRotation = maxRotation * rotationSize / slowRadius;
        }
        targetRotation *= rotation / rotationSize;


        float resultAngular = targetRotation - character->getAngularVelocity();
        resultAngular /= timeToTarget;

        float angularAcceleration = fabs(resultAngular);
        if (angularAcceleration > maxAngularAcceleration) {
            resultAngular /= angularAcceleration;
            resultAngular *= maxAngularAcceleration;
        }
        result->setAngular(resultAngular);
    }
    // todo implement the rest.

    
    result->setLinear(sf::Vector2f(0.0f, 0.0f));
    //result.angular = 0; // Assuming no angular acceleration for simplicity
    return result;
}


