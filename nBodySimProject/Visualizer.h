#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#pragma once
class Visualizer
{
public:
	//ctor takes the particle container, assigns it to the visualContainer and calls the visualizing
	Visualizer(const size_t i_numberOfParticles, const float i_scale, const float i_gravitationalConstant, const unsigned int i_edgeFreePixels, const unsigned int i_option, const float i_dt);
	void render();
	void synchSimAndVisualization();
private:
	std::vector<sf::CircleShape> m_circleContainer;
	void setUpCircleContainer();
	Simulation m_sim;
	const float m_scale{ 0.0F };
	sf::CircleShape m_tmpCircle;

	sf::VideoMode m_desktopMode{ sf::VideoMode::getDesktopMode() };
	const unsigned int m_edgeFreePixels{ 0U };
	const unsigned int m_option{ 0U };
	const float m_dt{ 0.0F };
};

