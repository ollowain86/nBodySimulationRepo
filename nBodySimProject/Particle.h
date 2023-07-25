#include "helperTypesAndFunctions.h"
#pragma once
class Particle
{
public:
	helpers::vec2f m_pos;
	helpers::vec2f m_vel;
	helpers::vec2f m_vel_half_dt;
	helpers::vec2f m_accel{0.0F, 0.0F};
	float m_mass;
private:
	
};

