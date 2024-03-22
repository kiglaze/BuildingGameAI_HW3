#ifndef SPRITECOLLECTION_H
#define SPRITECOLLECTION_H

#include <vector>
#include "Sprite.h"
#include "Crumb.h"

class SpriteCollection {
public:
    SpriteCollection(std::vector<Crumb>* crumbs);
    // Override the virtual destructor
    virtual ~SpriteCollection();

    void allSpritesDropCrumbs();
    void addSprite(Sprite* sprite);
    void removeSprite(Sprite* sprite);
    void drawAll(sf::RenderWindow& window);
    std::size_t getSpriteCount();
    void addStartingSprite(std::string textureFilePath, int isLeader);
    const std::vector<Sprite*>& getSprites() const;
    sf::Vector2f getCOMPosition();
    sf::Vector2f getCOMVelocityVector();
    void deleteMarkedSprites();

private:
    std::vector<Sprite*> sprites;
    std::vector<Crumb>* breadcrumbs;
    sf::Vector2f bc_position = sf::Vector2f(0.0f, 0.0f);

    float drop_timer = 15.f;
    int crumb_idx = 0;
};

#endif // SPRITECOLLECTION_H
