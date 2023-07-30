#include "Visualizer.h"
#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#include "Particle.h"
#include "TextHandler.h"
#include "RectangleHandler.h"
#include <string>

Visualizer::Visualizer(const size_t i_numberOfParticles, const float i_scale, const float i_gravitationalConstant, const unsigned int i_edgeFreePixels, const unsigned int i_option, const float i_dt) : m_sim(i_numberOfParticles, i_scale, i_gravitationalConstant, i_edgeFreePixels), m_scale(i_scale), m_edgeFreePixels(i_edgeFreePixels), m_option(i_option), m_dt(i_dt)
{
	//creates particle container with i_numberOfParticles, each particle has pos, vel and mass
	m_sim.setUpSelector(m_desktopMode.width, m_desktopMode.height, m_option);
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
	// Use the desktop resolution for the window size
	sf::RenderWindow window(sf::VideoMode(m_desktopMode.width, m_desktopMode.height), "N-Body Gravity Simulation");

	// Font for displaying the framerate
	sf::Font font;
	if (!font.loadFromFile("../nBodySimProject/arial.ttf"))
	{
		//do nothing
	}

	// Text object to display the framerate
	std::string tmpString = "FPS: ";
	unsigned int tmpCharSize = 20U;
	sf::Color tmpFillColor = sf::Color::White;
	sf::Color tmpOutlineColor = sf::Color::White;
	float tmpOutlineThickness = 0.0F;
	float tmpPosX = 0.0F;
	float tmpPosY = 0.0F;
	TextHandler textHandler;
	textHandler.setText(font, tmpString, tmpCharSize, tmpFillColor, tmpOutlineColor, tmpOutlineThickness, tmpPosX, tmpPosY);

	// Set clock for FPS calculation
	sf::Clock clock;
	sf::Time elapsed = sf::Time::Zero;
	unsigned int frameCount = 0;
	const sf::Time updateRate = sf::seconds(1.0f); // Update the FPS text every second

	window.setFramerateLimit(60);

	// THE SFML WHILE LOOPs
	while (window.isOpen()) {
		sf::Event event;
		// POLL EVENTS
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

		// CLEAR, DRAW, DISPLAY START
		window.clear(sf::Color::Black);

		for (size_t i = 0; i < m_circleContainer.size(); i++)
		{
			window.draw(m_circleContainer[i]);
		}

		window.draw(textHandler.getText());

		window.display();
		// CLEAR, DRAW, DISPLAY END

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

		m_sim.moveParticles(m_dt);
		synchSimAndVisualization();
		//m_sim.calcTotalEnergy();
		//m_sim.writeOutData();
	}
}