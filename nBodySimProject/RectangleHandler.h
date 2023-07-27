#include <SFML/Graphics.hpp>
#pragma once
class RectangleHandler
{
public:
	RectangleHandler();
	void setRectangle(const sf::Vector2f& i_pos, const sf::Color& i_outLineColor, const float outLineThickness, const sf::Color& i_fillColor, const sf::Vector2f& i_size);
	sf::RectangleShape& getRectangle();
private:
	sf::RectangleShape m_rectangle;
};

