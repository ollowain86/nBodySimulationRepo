#include "Visualizer.h"
#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#include "Particle.h"
#include "TextHandler.h"
#include "RectangleHandler.h"
#include <string>
#include <iostream>

Visualizer::Visualizer(const size_t i_numberOfParticles, const double i_scale, const double i_gravitationalConstant, const unsigned int i_edgeFreePixels, const unsigned int i_option, const double i_dt) : m_sim(i_numberOfParticles, i_scale, i_gravitationalConstant, i_edgeFreePixels, i_option), m_scale(i_scale), m_edgeFreePixels(i_edgeFreePixels), m_option(i_option), m_dt(i_dt)
{
	//creates particle container with i_numberOfParticles, each particle has pos, vel and mass
	m_sim.setUpSelector(m_desktopMode.width, m_desktopMode.height);
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

		m_tmpCircle.setRadius(5.0);
		m_tmpCircle.setPosition((tmpParticle.m_pos.x - m_tmpCircle.getRadius()/m_scale), (tmpParticle.m_pos.y - m_tmpCircle.getRadius() / m_scale));
		m_tmpCircle.setFillColor(sf::Color::White);
		m_tmpCircle.setOutlineColor(sf::Color::Blue);
		m_tmpCircle.setOutlineThickness(2.0);

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
	double tmpOutlineThickness = 0.0;
	double tmpPosX = 0.0;
	double tmpPosY = 0.0;
	TextHandler textHandler;
	textHandler.setText(font, tmpString, tmpCharSize, tmpFillColor, tmpOutlineColor, tmpOutlineThickness, tmpPosX, tmpPosY);

	// Set clock for FPS calculation
	sf::Clock clock;
	sf::Time elapsed = sf::Time::Zero;
	unsigned int frameCount = 0;
	const sf::Time updateRate = sf::seconds(1.0); // Update the FPS text every second

	//window.setFramerateLimit(60);

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
			// Check for mouse scroll events
			if (event.type == sf::Event::MouseWheelScrolled)
			{
				if (event.mouseWheelScroll.delta < 0)
				{
					m_zoomFactor *= 1.1; // Increase zoom on scroll up
				}
				else if (event.mouseWheelScroll.delta > 0)
				{
					m_zoomFactor *= 0.9; // Decrease zoom on scroll down
				}
				else
				{
					//do nothing
				}
				// Clamp the zoom factor to reasonable values
				m_zoomFactor = std::max(0.1, std::min(100.0, m_zoomFactor));
				// 
				//textHandler.setCharacterSize(static_cast<unsigned int>(textHandler.getText().getCharacterSize() * m_zoomFactor));
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

		// Set the view
		sf::View view = window.getDefaultView();
		view.zoom(m_zoomFactor);
		window.setView(view);
		// CLEAR, DRAW, DISPLAY END

		// Calculate framerate
		elapsed += clock.restart();
		frameCount++;
		if (elapsed >= updateRate)
		{
			double fps = static_cast<double>(frameCount) / elapsed.asSeconds();
			textHandler.setString("FPS: " + std::to_string(static_cast<int>(fps)));

			elapsed = sf::Time::Zero;
			frameCount = 0;
		}

		m_sim.moveParticles(m_dt);
		synchSimAndVisualization();
		m_sim.writeOutData();
	}
}