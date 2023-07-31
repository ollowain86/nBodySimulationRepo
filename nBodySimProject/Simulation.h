#include <vector>
#include "Particle.h"
#include "helperTypesAndFunctions.h"

#pragma once
class Simulation
{
public:
	Simulation(const size_t i_numOfParticles, const float i_scale, const float i_gravitationalConstant, const unsigned int i_edgeFreePixels);
	//return the m_particleContainer
	const std::vector<Particle>& getParticleContainer() const;
	//moves the particles for dt
	void moveParticles(const float i_dt);
	void writeOutData();
	// TESTING AND DEBUGGING
	// i_option: 0U = normal, 1U predefined 2 particle system, 2U predefined 3 particle system
	void setUpSelector(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr, const unsigned int i_option);
	void setUpTwoParticle();
	void setUpThreeParticle();
	bool floatEqual(const float a, const float b);
	// calculate total potential Energy
	float calcTotalPotentialEnergy();
	// calculate total kinetic Energy
	float calcTotalKineticEnergy();
private:
	//######## METHODS ########
	void leapfrogUpdate(const float i_dt);
	void calculateAcceleration(Particle& particle);
	//calculates orbital velocity with v_c = sqrt(G*M/r);
	void calcOrbitalSpeed(Particle& particle, const float i_massWithinR_i);
	// calcDistance between two particles
	float calcDistance(const Particle& particleA, const Particle& particleB);
	// calculate length of vec2f
	float calcLength(const helpers::vec2f i_2dVec);
	//calculates M_i (mass inside r_i)
	float massWithinRadiusCalculator(const Particle& particle);
	// particles positions are set randomly into a circular shape with their radius (no vel, no accel assigned here), but m_particleContainer gets its initial size
	void setUpCircularShape(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr, const float i_maxRadius);
	// assigns pos, vel and mass to the particles
	void setUpSimulation(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr);
	//######## PARAMETERS ########
	// number of Particles in the Sim
	const size_t m_numberOfParticles{ 0U };
	const float m_scale{ 0.0F };
	// particle container
	std::vector<Particle> m_particleContainer;
	// simulation parameters
	float m_deltaX{ 0.0F };
	float m_deltaY{ 0.0F };
	float m_distance{ 0.0F };
	float m_distanceSqrd{ 0.0F };
	float m_plummerRadius{ 1.0F };
	const float m_gravitationalConstant{ 1.0F };
	float m_accelMagnitude{ 0.0F };
	float m_accelerationX{ 0.0F };
	float m_accelerationY{ 0.0F };
	//important for visualizing and particle distribution
	const unsigned int m_edgeFreePixels{ 0U };
};

