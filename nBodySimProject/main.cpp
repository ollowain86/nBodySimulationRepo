#include <iostream>
#include "SFML/Graphics.hpp"
#include "Simulation.h"
#include <vector>
#include "Particle.h"
#include "SFML/Graphics.hpp"


int main()
{
    //number of particles
    size_t numOfParticles{ 2U };
    Simulation simulation(numOfParticles);


    std::vector<Particle> particleContainer;

    Particle tmpParticle;
    for (size_t i = 0; i < numOfParticles; i++)
    {
        tmpParticle.setMass(10.0F);
        tmpParticle.setPos({(i+1)*100.0F, (i+1) *100.0F});
        tmpParticle.setVel({ (i + 1) *3.0F, (i + 1) *3.0F });

        particleContainer.push_back(tmpParticle);
    }

    // SFML
    sf::RenderWindow window(sf::VideoMode(600, 600), "SFML Circle");
    std::vector<sf::CircleShape> circleContainer;
    sf::CircleShape tmpCircle;
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::Black);

        for (size_t i = 0; i < particleContainer.size(); i++)
        {
            tmpCircle.setFillColor(sf::Color::White);
            tmpCircle.setOutlineColor(sf::Color::Blue);
            tmpCircle.setOutlineThickness(2.0F);
            tmpCircle.setRadius(10.0F);
            tmpCircle.setPosition(particleContainer[i].getPos().x- tmpCircle.getRadius(), particleContainer[i].getPos().y- tmpCircle.getRadius());
            circleContainer.push_back(tmpCircle);
        }

        for (size_t i = 0; i < circleContainer.size(); i++)
        {
            window.draw(circleContainer[i]);
        }

        window.display();
    }

    return 0;
}