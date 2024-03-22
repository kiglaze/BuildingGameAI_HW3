#ifndef LOOK_WHERE_YOURE_GOING_H
#define LOOK_WHERE_YOURE_GOING_H

#include "Align.h" // Include the base class header
#include "Kinematic.h"
#include "Sprite.h"
#include "SteeringData.h"

class LookWhereYoureGoing : public Align {
public:
    // Constructor with initialization list, reusing Align's constructor
    LookWhereYoureGoing(Kinematic* target, Sprite* character);

    // Override the virtual destructor
    virtual ~LookWhereYoureGoing() {}

    // Override execute() method from SteeringBehavior to include facing logic
    void execute(float timeDelta) override;

    // Override calculateAcceleration() method from SteeringBehavior to calculate the facing steering
    SteeringData calculateAcceleration() override;

    // Override calculateAccelerationPointer() method from SteeringBehavior for pointer version
    SteeringData* calculateAccelerationPointer() override;

protected:

private:
    // Private methods or member variables specific to Face, if needed
};

#endif // LOOK_WHERE_YOURE_GOING_H
