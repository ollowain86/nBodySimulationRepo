#include "helperTypesAndFunctions.h"
#pragma once
class Particle
{
public:
	void setMass(const float i_mass);
	void setPos(const helpers::vec2f& i_Pos);
	void setVel(const helpers::vec2f& i_Vel);
	float getMass();
	helpers::vec2f getPos();
	helpers::vec2f getVel();
private:
	helpers::vec2f m_pos;
	helpers::vec2f m_vel;
	float m_mass;
};

