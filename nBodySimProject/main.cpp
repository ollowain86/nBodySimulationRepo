#include <iostream>
#include "SFML/Graphics.hpp"

int main()
{
    sf::RenderWindow window(sf::VideoMode(800, 600), "SFML Circle");

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

        // Create a circle shape
        sf::CircleShape circle(50.0f);
        circle.setFillColor(sf::Color::Red);
        circle.setPosition(375.0f, 275.0f);

        window.draw(circle);

        window.display();
    }
    return 0;
}