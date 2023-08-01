#include "TextHandler.h"
#include <SFML/Graphics.hpp>

TextHandler::TextHandler()
{
}

void TextHandler::setText(const sf::Font& i_font, const std::string& i_string, const unsigned int i_charSize, const sf::Color& i_fillColor, const sf::Color& i_outLineColor, const double i_outLineThickness, const double posX, const double posY)
{
    m_text.setFont(i_font);
    m_text.setString(i_string);
    m_text.setCharacterSize(i_charSize);
    m_text.setFillColor(i_fillColor);
    m_text.setOutlineColor(i_outLineColor); // Outline color
    m_text.setOutlineThickness(i_outLineThickness);
    m_text.setPosition(posX, posY);
}

void TextHandler::setString(const std::string& i_string)
{
    m_text.setString(i_string);
}

sf::Text& TextHandler::getText()
{
    return m_text;
}

// Setter for position
void TextHandler::setPosition(const double posX, const double posY)
{
    m_text.setPosition(posX, posY);
}

// Setter for character size
void TextHandler::setCharacterSize(const unsigned int i_charSize)
{
    m_text.setCharacterSize(i_charSize);
}