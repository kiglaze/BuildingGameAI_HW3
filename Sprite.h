#ifndef SPRITE_H
#define SPRITE_H

#include <SFML/Graphics.hpp>
#include <string>
#include <vector>
#include "Kinematic.h"
#include "Crumb.h"

class SpriteCollection; // Forward declaration if SpriteCollection is used

class Sprite : public Kinematic {
public:
    Sprite(const std::string& textureFile, float posX, float posY, float orient, sf::Vector2f vel, float rot, int isLeader, std::vector<Crumb>* crumbs);

    void dropSomeCrumbs();

    sf::Vector2f getVelocityVector() const override;
    float getSpriteWidth();
    float getSpriteHeight();

    void move(float dx, float dy);
    void moveRight(float diff);
    void moveDown(float diff);
    void moveLeft(float diff);
    void moveUp(float diff);

    void setScale(float scaleX, float scaleY);
    void setVelocityVector(float velX, float velY) override;
    void draw(sf::RenderWindow& window);
    sf::Vector2f getPosition() const override;
    void setPosition(float posX, float posY) override;


    float mapToRange(float rotation);
    void setDirection(float newDirection) override;
    void turnRight();
    float getDirection() const override;

    void rotate(float x);
    void setRotation(float x);
    float getRotation();
    float getAngularVelocity() const override;
    void setAngularVelocity(float w) override;

    int getHasStarted();
    void setHasStarted(int hasStartedVal);

    int getHasReachedCorner();
    int getIsLeader();
    void markForDeletion();
    int shouldBeDeleted();

    void moveWithVelocityVector(float timeDelta);
    void moveAccordingToDirection(float timeDelta, float screenWidth, float screenHeight);
    void rotateIfNeeded(float screenWidth, float screenHeight);

private:
    sf::Texture texture;
    sf::Sprite sprite;
    //float direction;
    //float angular_velocity;
    int has_started;
    int has_reached_corner;
    int is_leader;
    int should_delete;
    //sf::Vector2f velocity;

    std::vector<Crumb>* breadcrumbs;
    float drop_timer = 15.f;
    int crumb_idx = 0;
    sf::Vector2f bc_position = sf::Vector2f(0.0f, 0.0f);
};

#endif // SPRITE_H
