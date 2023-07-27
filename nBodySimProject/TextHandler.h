#include <SFML/Graphics.hpp>
#include <string>
#pragma once
class TextHandler
{
public:
	TextHandler();
	void setText(const sf::Font& i_font, const std::string& i_string, const unsigned int i_charSize, const sf::Color& i_fillColor, const sf::Color& i_outLineColor, const float i_outLineThickness, const float posX, const float posY);
	void setString(const std::string& i_string);
	sf::Text& getText();
private:
	sf::Text m_text;
};

