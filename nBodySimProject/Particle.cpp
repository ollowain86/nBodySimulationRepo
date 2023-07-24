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

float Particle::getMass()
{
	return m_mass;
}

helpers::vec2f Particle::getVel()
{
	return m_vel;
}
helpers::vec2f Particle::getPos()
{
	return m_pos;
}