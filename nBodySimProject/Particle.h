#include "helperTypesAndFunctions.h"
#pragma once
class Particle
{
public:
	helpers::vec2d_double m_pos{ 0.0, 0.0 };
	helpers::vec2d_double m_vel{ 0.0, 0.0 };
	helpers::vec2d_double m_vel_half_dt{ 0.0, 0.0 };
	helpers::vec2d_double m_accel{0.0, 0.0};
	double m_mass;
	double m_velScalar;
	// radius around the chosen origin
	double m_radius{ 0.0 };
private:
	
};

