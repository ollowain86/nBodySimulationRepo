#include <vector>
#include "Particle.h"
#include "helperTypesAndFunctions.h"

#pragma once
class Simulation
{
public:
	Simulation(const size_t i_numOfParticles, const double i_scale, const double i_gravitationalConstant, const unsigned int i_edgeFreePixels);
	//return the m_particleContainer
	const std::vector<Particle>& getParticleContainer() const;
	//moves the particles for dt
	void moveParticles(const double i_dt);
	void writeOutData();
	// TESTING AND DEBUGGING
	// i_option: 0U = normal, 1U predefined 2 particle system, 2U predefined 3 particle system
	void setUpSelector(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr, const unsigned int i_option);
	void setUpTwoParticle();
	double calcTotalEnergy();
	void setUpThreeParticle();
	bool doubleEqual(const double a, const double b);
	// calculate total potential Energy
	double calcTotalPotentialEnergy();
	// calculate total kinetic Energy
	double calcTotalKineticEnergy();
private:
	//######## METHODS ########
	void leapfrogUpdate(const double i_dt);
	void calculateAcceleration(Particle& particle);
	//calculates orbital velocity with v_c = sqrt(G*M/r);
	void calcOrbitalSpeed(Particle& particle, const double i_massWithinR_i);
	// calcDistance between two particles
	double calcDistance(const Particle& particleA, const Particle& particleB);
	// calculate length of vec2f
	double calcLength(const helpers::vec2d_double i_2dVec);
	//calculates M_i (mass inside r_i)
	double massWithinRadiusCalculator(const Particle& particle);
	// particles positions are set randomly into a circular shape with their radius (no vel, no accel assigned here), but m_particleContainer gets its initial size
	void setUpCircularShape(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr, const double i_maxRadius);
	// assigns pos, vel and mass to the particles
	void setUpSimulation(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr);
	//######## PARAMETERS ########
	// number of Particles in the Sim
	const size_t m_numberOfParticles{ 0U };
	const double m_scale{ 0.0 };
	// particle container
	std::vector<Particle> m_particleContainer;
	// simulation parameters
	double m_deltaX{ 0.0 };
	double m_deltaY{ 0.0 };
	double m_distance{ 0.0 };
	double m_distanceSqrd{ 0.0 };
	double m_plummerRadius{ 1.0 };
	const double m_gravitationalConstant{ 1.0 };
	double m_accelMagnitude{ 0.0 };
	double m_accelerationX{ 0.0 };
	double m_accelerationY{ 0.0 };
	//important for visualizing and particle distribution
	const unsigned int m_edgeFreePixels{ 0U };
};

