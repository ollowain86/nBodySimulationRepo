#include "Visualizer.h"
#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#include "Particle.h"
#include "TextHandler.h"
#include "RectangleHandler.h"
#include <string>

Visualizer::Visualizer(const size_t i_numberOfParticles, const float i_scale) : m_sim(i_numberOfParticles, i_scale), m_scale(i_scale)
{
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

		m_tmpCircle.setRadius(5.0F);
		m_tmpCircle.setPosition((tmpParticle.m_pos.x - m_tmpCircle.getRadius()/m_scale), (tmpParticle.m_pos.y - m_tmpCircle.getRadius() / m_scale));
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
		m_circleContainer[i].setPosition((tmpParticle.m_pos.x - tmpCircle.getRadius()) / m_scale, (tmpParticle.m_pos.y - tmpCircle.getRadius()) / m_scale);
	}
}

void Visualizer::render()
{
	// SFML
	sf::RenderWindow window(sf::VideoMode(1900, 1000), "N-Body Gravity Simulation");

	// Font for displaying the framerate
	sf::Font font;
	if (!font.loadFromFile("../nBodySimProject/arial.ttf"))
	{
		// Handle font loading error
		// Replace "arial.ttf" with the path to your desired font file
	}

	// Text object to display the framerate
	TextHandler textHandler;
	textHandler.setText(font, "FPS: ", 20U, sf::Color::White, sf::Color::White, 0.0F, 10.0F, 10.0F);

	sf::Clock clock;
	sf::Time elapsed = sf::Time::Zero;
	unsigned int frameCount = 0;
	const sf::Time updateRate = sf::seconds(1.0f); // Update the FPS text every second

	window.setFramerateLimit(60);

	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			if (event.type == sf::Event::Closed)
			{
				window.close();
			}
			// Check for the Escape key press
			if (event.type == sf::Event::KeyPressed)
			{
				if (event.key.code == sf::Keyboard::Escape)
				{
					window.close();
				}
			}
		}

		window.clear(sf::Color::Black);

		for (size_t i = 0; i < m_circleContainer.size(); i++)
		{
			window.draw(m_circleContainer[i]);
		}

		window.draw(textHandler.getText());

		window.display();

		// Calculate framerate
		elapsed += clock.restart();
		frameCount++;
		if (elapsed >= updateRate)
		{
			float fps = static_cast<float>(frameCount) / elapsed.asSeconds();
			textHandler.setString("FPS: " + std::to_string(static_cast<int>(fps)));

			elapsed = sf::Time::Zero;
			frameCount = 0;
		}

		m_sim.moveParticles(1.0F);
		synchSimAndVisualization();
		//m_sim.calcTotalEnergy();
		//m_sim.writeOutData();
	}
}