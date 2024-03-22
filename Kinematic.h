// Kinematic.h
#ifndef KINEMATIC_H
#define KINEMATIC_H

#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include "SteeringData.h"

class Kinematic {
public:
    Kinematic(sf::Vector2f pos, float orient, sf::Vector2f vel, float rot);
    virtual ~Kinematic() {}

    virtual sf::Vector2f getPosition() const;
    virtual void setPosition(float posX, float posY);
    virtual float getDirection() const;
    float mapToRange(float rotation);
    virtual void setDirection(float newDirection);
    virtual sf::Vector2f getVelocityVector() const;
    virtual void setVelocityVector(float valX, float valY);
    virtual float getAngularVelocity() const;
    virtual void setAngularVelocity(float w);
    void update(sf::Vector2f positionVal, float orientationVal, sf::Vector2f velocityVal, float rotationVal);
    void update(SteeringData sd, float deltaTime);
protected:
    sf::Vector2f position;
    float direction; // orientation
    sf::Vector2f velocity;
    float angular_velocity; // rotation

    float maxWindowX = 640.0;
    float maxWindowY = 480.0;
};

#endif // KINEMATIC_H
