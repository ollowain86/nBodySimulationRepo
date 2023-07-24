#include "Particle.h"
#include "helperTypesAndFunctions.h"

void Particle::setMass(const float i_mass)
{
	m_mass = i_mass;
}

void Particle::setVel(const helpers::vec2f& i_Vel)
{
	m_vel = i_Vel;
}

void Particle::setPos(const helpers::vec2f& i_Pos)
{
	m_pos = i_Pos;
}