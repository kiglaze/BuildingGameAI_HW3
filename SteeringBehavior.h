#ifndef STEERING_BEHAVIOR_H
#define STEERING_BEHAVIOR_H

#include "SteeringData.h"
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>

class SteeringBehavior {

public:
    // Virtual destructor to ensure derived class destructors are called
    virtual ~SteeringBehavior() {}

    float getLengthOfVector(sf::Vector2f vectorInput);

    void normalizeVector(sf::Vector2f& vectorInput);

    // Pure virtual function for executing the behavior
    virtual void execute(float timeDelta) = 0;

    virtual SteeringData calculateAcceleration() = 0;
    virtual SteeringData* calculateAccelerationPointer() = 0;
};

#endif // STEERING_BEHAVIOR_H
