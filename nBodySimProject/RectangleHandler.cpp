#include "RectangleHandler.h"

RectangleHandler::RectangleHandler()
{

}

void RectangleHandler::setRectangle(const sf::Vector2f& i_pos, const sf::Color& i_outLineColor, const float outLineThickness, const sf::Color& i_fillColor, const sf::Vector2f& i_size)
{
    m_rectangle.setPosition(i_pos);
    m_rectangle.setOutlineColor(i_outLineColor);
    m_rectangle.setOutlineThickness(outLineThickness);
    m_rectangle.setFillColor(i_fillColor);
    m_rectangle.setSize(i_size);
}

sf::RectangleShape& RectangleHandler::getRectangle()
{
    return m_rectangle;
}