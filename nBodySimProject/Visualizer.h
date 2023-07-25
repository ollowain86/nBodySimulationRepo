#include "SFML/Graphics.hpp"
#include <vector>
#include "Simulation.h"
#pragma once
class Visualizer
{
public:
	//ctor takes the particle container, assigns it to the visualContainer and calls the visualizing
	Visualizer(const size_t i_numberOfParticles);
	void render();
	void synchSimAndVisualization();
private:
	std::vector<sf::CircleShape> m_circleContainer;
	void setUpCircleContainer();
	Simulation m_sim;
	sf::CircleShape m_tmpCircle;
};

