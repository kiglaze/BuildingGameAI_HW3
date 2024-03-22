#ifndef ALIGN_H
#define ALIGN_H

#include "SteeringData.h"
#include <SFML/Graphics.hpp>
#include "SteeringBehavior.h"
#include "Sprite.h" // Make sure this exists and is correctly implemented

class Align : public SteeringBehavior {
public:
    // Constructor with initialization list
    Align(Kinematic* target, Sprite* character);

    // Override the virtual destructor
    virtual ~Align();

    void setTargetRadius(float targetRad);

    void setMaxAngularAcceleration(float maxAngAcc);
    void setMaxRotation(float maxRot);
    void setSlowRadius(float slowRadiusVal);
    void setTimeToTarget(float timeToTargetVal);

    float mapToRange(float rotation);

    // Override execute() method from SteeringBehavior
    void execute(float timeDelta) override;

    // Override calculateAcceleration() method from SteeringBehavior
    SteeringData calculateAcceleration() override;

    // Override calculateAcceleration() method from SteeringBehavior
    SteeringData* calculateAccelerationPointer() override;

protected:
    Kinematic* target;
    Sprite* character;
    float targetRadius = 15;

    float maxAngularAcceleration = .01;
    float maxRotation = 0.1;
    float slowRadius = 45;
    float timeToTarget = 6;
};

#endif // ALIGN_H
