#ifndef ARRIVE_H
#define ARRIVE_H

#include "SteeringData.h"
#include <SFML/Graphics.hpp>
#include "SteeringBehavior.h"
#include "Sprite.h" // Make sure this exists and is correctly implemented
#include "Kinematic.h"

class Arrive : public SteeringBehavior {
public:
    // Constructor with initialization list
    Arrive(Kinematic* target, Sprite* character);

    // Override the virtual destructor
    virtual ~Arrive();

    void setTargetRadius(float targetRadiusVal);
    void setSlowRadius(float slowRadiusVal);

    // Override execute() method from SteeringBehavior
    void execute(float timeDelta) override;

    // Override calculateAcceleration() method from SteeringBehavior
    SteeringData calculateAcceleration() override;

    // Override calculateAcceleration() method from SteeringBehavior
    SteeringData* calculateAccelerationPointer() override;
protected:
    Kinematic* target;
    Sprite* character;
private:
    // Private member variables
    float maxSpeed = 0.1;
    float targetRadius = 50;
    float slowRadius = 60;
    float maxAcceleration = 0.002;
    float timeToTarget = 1.0;
};

#endif // ARRIVE_H
