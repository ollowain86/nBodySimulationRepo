#include <vector>
#include "Particle.h"
#pragma once
class Simulation
{
public:
	Simulation(const size_t i_numOfParticles, const float i_scale);
	// assigns pos, vel and mass to the particles
	void setUpSimulation();
	//return the m_particleContainer
	const std::vector<Particle>& getParticleContainer() const;
	//moves the particles for dt
	void moveParticles(const float i_dt);
	void calcTotalEnergy();
	void writeOutData();
private:
	//######## METHODS ########
	void leapfrogUpdate(const float i_dt);
	void calculateAcceleration(Particle& particle);
	//######## PARAMETERS ########
	// number of Particles in the Sim
	const size_t m_numberOfParticles{ 0U };
	const float m_scale{ 0.0F };
	// particle container
	std::vector<Particle> m_particleContainer;
};

