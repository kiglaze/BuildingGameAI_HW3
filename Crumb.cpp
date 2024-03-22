#include <SFML/Graphics.hpp>
#include "Crumb.h"
//Breadcrumb class

Crumb::Crumb(int id)
{
    //set initial position and size breadcrumbs   
    this->id = id;         
    this->setRadius(10.f);
    this->setFillColor(sf::Color(0, 0, 255, 255));
    this->setPosition(-100, -100);
}

//tell breadcrumb to render self, using current render window
void Crumb::draw(sf::RenderWindow* window)
{
    window->draw(*this);
}

//set position of breadcrumb
void Crumb::drop(float x, float y)
{
    this->setPosition(x, y);
}

//set position of breadcrumb
void Crumb::drop(sf::Vector2f position)
{
    this->setPosition(position);
}

