#include "helperTypesAndFunctions.h"
#pragma once
class Particle
{
public:
	helpers::vec2f m_pos{ 0.0F, 0.0F };
	helpers::vec2f m_vel{ 0.0F, 0.0F };
	helpers::vec2f m_vel_half_dt{ 0.0F, 0.0F };
	helpers::vec2f m_accel{0.0F, 0.0F};
	float m_mass;
private:
	
};

