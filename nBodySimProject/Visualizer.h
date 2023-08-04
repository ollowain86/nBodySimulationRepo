#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#pragma once
class Visualizer
{
public:
	//ctor takes the particle container, assigns it to the visualContainer and calls the visualizing
	Visualizer(const size_t i_numberOfParticles, const double i_scale, const double i_gravitationalConstant, const unsigned int i_edgeFreePixels, const unsigned int i_option, const double i_dt);
	void render();
	void synchSimAndVisualization();
private:
	std::vector<sf::CircleShape> m_circleContainer;
	void setUpCircleContainer();
	Simulation m_sim;
	const double m_scale{ 0.0 };
	sf::CircleShape m_tmpCircle;

	sf::VideoMode m_desktopMode{ sf::VideoMode::getDesktopMode() };
	const unsigned int m_edgeFreePixels{ 0U };
	const unsigned int m_option{ 0U };
	const double m_dt{ 0.0 };
	double m_zoomFactor{ 1.0 }; //for visual zoom
};

