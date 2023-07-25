#include "Visualizer.h"
#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#include "Particle.h"

Visualizer::Visualizer(const size_t i_numberOfParticles)
{
	m_sim.setNumOfParticles(i_numberOfParticles);
	//creates particle container with i_numberOfParticles, each particle has pos, vel and mass
	m_sim.setUpSimulation();
	// takes the particle from the simulation and fills the m_circleContainer
	setUpCircleContainer();
}

void Visualizer::setUpCircleContainer()
{
	//for readability
	Particle tmpParticle;

	for (size_t i = 0; i < m_sim.getParticleContainer().size(); i++)
	{
		tmpParticle = m_sim.getParticleContainer().at(i);

		m_tmpCircle.setRadius(10.0F);
		m_tmpCircle.setPosition(tmpParticle.m_pos.x - m_tmpCircle.getRadius(), tmpParticle.m_pos.y - m_tmpCircle.getRadius());
		m_tmpCircle.setFillColor(sf::Color::White);
		m_tmpCircle.setOutlineColor(sf::Color::Blue);
		m_tmpCircle.setOutlineThickness(2.0F);

		m_circleContainer.push_back(m_tmpCircle);
	}
}

void Visualizer::synchSimAndVisualization()
{
	sf::CircleShape tmpCircle;
	Particle tmpParticle;
	for (size_t i = 0; i < m_sim.getParticleContainer().size(); i++)
	{
		tmpParticle = m_sim.getParticleContainer().at(i);
		m_circleContainer[i].setPosition(tmpParticle.m_pos.x - tmpCircle.getRadius(), tmpParticle.m_pos.y - tmpCircle.getRadius());
	}
}

void Visualizer::render()
{
	// SFML
	sf::RenderWindow window(sf::VideoMode(1900, 1000), "N-Body Gravity Simulation");

	window.setFramerateLimit(60);

	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed)
				window.close();
		}

		window.clear(sf::Color::Black);

		for (size_t i = 0; i < m_circleContainer.size(); i++)
		{
			window.draw(m_circleContainer[i]);
		}

		window.display();

		m_sim.moveParticles(0.1F);
		synchSimAndVisualization();
		m_sim.calcTotalEnergy();
	}
}