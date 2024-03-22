#ifndef FACE_H
#define FACE_H

#include "Align.h" // Include the base class header
#include "Kinematic.h"
#include "Sprite.h"
#include "SteeringData.h"

class Face : public Align {
public:
    // Constructor with initialization list, reusing Align's constructor
    Face(Kinematic* target, Sprite* character);

    // Override the virtual destructor
    virtual ~Face() {}

    // Override execute() method from SteeringBehavior to include facing logic
    void execute(float timeDelta) override;

    // Override calculateAcceleration() method from SteeringBehavior to calculate the facing steering
    SteeringData calculateAcceleration() override;

    // Override calculateAccelerationPointer() method from SteeringBehavior for pointer version
    SteeringData* calculateAccelerationPointer() override;

protected:
    Kinematic* target = nullptr;
    // Additional methods or member variables specific to Face, if needed

private:
    // Private methods or member variables specific to Face, if needed
};

#endif // FACE_H
