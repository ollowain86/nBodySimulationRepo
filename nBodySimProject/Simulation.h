#include <vector>
#include "Particle.h"
#pragma once
class Simulation
{
public:
	// set number of Particles
	void setNumOfParticles(const size_t i_numOfParticles);
	// assigns pos, vel and mass to the particles
	void setUpSimulation();
	//return the m_particleContainer
	const std::vector<Particle>& getParticleContainer() const;
	//moves the particles for dt
	void moveParticles(const float i_dt);
	void calcTotalEnergy();
private:
	//######## METHODS ########
	void leapfrogUpdate(const float i_dt);
	void calculateAcceleration(Particle& particle);
	//######## PARAMETERS ########
	// number of Particles in the Sim
	size_t m_numberOfParticles{ 0U };
	// particle container
	std::vector<Particle> m_particleContainer;
};

