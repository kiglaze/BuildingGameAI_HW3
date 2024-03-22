#ifndef CRUMB_HPP
#define CRUMB_HPP

#include <SFML/Graphics.hpp>

class Crumb : public sf::CircleShape {
public:
    explicit Crumb(int id);

    void draw(sf::RenderWindow* window);
    void drop(float x, float y);
    void drop(sf::Vector2f position);

private:
    int id;
};

#endif // CRUMB_HPP
