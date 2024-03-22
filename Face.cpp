#include "Face.h"
#include "Kinematic.h"
#include "Sprite.h"
#include "SteeringData.h"
#include <iostream>

// Constructor with initialization list, reusing Align's constructor
Face::Face(Kinematic* target, Sprite* character) : Align(target, character) {}

// Override execute() method from SteeringBehavior to include facing logic
void Face::execute(float timeDelta) {
    Align::execute(timeDelta);
}

// Override calculateAcceleration() method from SteeringBehavior to calculate the facing steering
SteeringData Face::calculateAcceleration() {
    return Align::calculateAcceleration();
}

// Override calculateAccelerationPointer() method from SteeringBehavior for pointer version
SteeringData* Face::calculateAccelerationPointer() {
    SteeringData* result = new SteeringData();
    sf::Vector2f parentTargetPosVect = Align::target->getPosition();
    sf::Vector2f direction = parentTargetPosVect - character->getPosition();
    if (getLengthOfVector(direction) == 0) {
        delete result;
        result = nullptr;
    } else {
        float targetOrientation = atan2(direction.y, direction.x) * (180.0 / M_PI);
        target = new Kinematic(parentTargetPosVect, targetOrientation, sf::Vector2f(0.0f, 0.0f), 0.0f);
        
        Align alignBehavior(target, character);
        alignBehavior.setTargetRadius(this->targetRadius);
        result = alignBehavior.calculateAccelerationPointer();
    }
    return result;
}
