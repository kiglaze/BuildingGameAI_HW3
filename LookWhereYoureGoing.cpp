#include "Face.h"
#include "Kinematic.h"
#include "Sprite.h"
#include "SteeringData.h"
#include <iostream>
#include "LookWhereYoureGoing.h"

// Constructor with initialization list, reusing Align's constructor
LookWhereYoureGoing::LookWhereYoureGoing(Kinematic* target, Sprite* character) : Align(target, character) {}

// Override execute() method from SteeringBehavior to include facing logic
void LookWhereYoureGoing::execute(float timeDelta) {
    Align::execute(timeDelta);
}

// Override calculateAcceleration() method from SteeringBehavior to calculate the facing steering
SteeringData LookWhereYoureGoing::calculateAcceleration() {
    return Align::calculateAcceleration();
}

// Override calculateAccelerationPointer() method from SteeringBehavior for pointer version
SteeringData* LookWhereYoureGoing::calculateAccelerationPointer() {
    SteeringData* result = new SteeringData();
    sf::Vector2f characterVelVect = character->getVelocityVector();
    if (getLengthOfVector(characterVelVect) == 0) {
        delete result;
        result = nullptr;
    } else {
        float targetOrientation = atan2(characterVelVect.y, characterVelVect.x) * (180.0 / M_PI);
        target->setDirection(targetOrientation);
        Align alignBehavior(target, character);
        alignBehavior.setTargetRadius(this->targetRadius);
        alignBehavior.setSlowRadius(5);
        alignBehavior.setTimeToTarget(1);
        alignBehavior.setMaxAngularAcceleration(0.05);
        result = alignBehavior.calculateAccelerationPointer();
    }
    return result;
}
