#pragma once
class Simulation
{
public:
	Simulation();
	Simulation(const size_t i_numberOfParticles);
private:
	size_t m_numberOfParticles{ 0U };
};

