#define _USE_MATH_DEFINES //for M_PI
#include <cmath>
#include "Simulation.h"
#include "Particle.h"
#include <vector>
#include <iostream>

void Simulation::setNumOfParticles(const size_t i_numOfParticles)
{
    m_numberOfParticles = i_numOfParticles;
}

// assigns pos, vel and mass to the particles
void Simulation::setUpSimulation()
{
    Particle tmpParticle;
    for (size_t i = 0; i < m_numberOfParticles; i++)
    {
        if (i==0)
        {
            tmpParticle.m_mass = 1000.0F;
            tmpParticle.m_pos.x = 2.0F;
            tmpParticle.m_pos.y = 4.0F;
            tmpParticle.m_vel.x = 3.0F;
            tmpParticle.m_vel.y = -1.0F;
            tmpParticle.m_accel.x = 0.0F;
            tmpParticle.m_accel.y = 0.0F;
        }
        if (i == 1)
        {
            tmpParticle.m_mass = 2000.0F;
            tmpParticle.m_pos.x = 6.0F;
            tmpParticle.m_pos.y = 1.0F;
            tmpParticle.m_vel.x = 0.0F;
            tmpParticle.m_vel.y = 1.0F;
            tmpParticle.m_accel.x = 0.0F;
            tmpParticle.m_accel.y = 0.0F;
        }

        if (i == 2)
        {
            tmpParticle.m_mass = 3000.0F;
            tmpParticle.m_pos.x = 1.0F;
            tmpParticle.m_pos.y = 1.0F;
            tmpParticle.m_vel.x = 1.0F;
            tmpParticle.m_vel.y = 0.0F;
            tmpParticle.m_accel.x = 0.0F;
            tmpParticle.m_accel.y = 0.0F;
        }

        m_particleContainer.push_back(tmpParticle);
    }

    // calc accel initially
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        calculateAcceleration(m_particleContainer[i]);
    }
}

const std::vector<Particle>& Simulation::getParticleContainer() const
{
    return m_particleContainer;
}

void Simulation::calculateAcceleration(Particle& particle)
{
    float deltaX{ 0.0F };
    float deltaY{ 0.0F };
    float distanceSquared{ 0.0F };
    float alpha{ 0.0F };
    const float gravitationalConstant{ 1.0F };
    for (const Particle& otherParticle : m_particleContainer)
    {
        if (&particle != &otherParticle)  // Ensure not calculating for the same particle
        {
            // Calculate distance between the two particles
            deltaX = otherParticle.m_pos.x - particle.m_pos.x;
            deltaY = otherParticle.m_pos.y - particle.m_pos.y;
            distanceSquared = deltaX * deltaX + deltaY * deltaY;

            // Calculate the gravitational force magnitude
            float forceMagnitude = gravitationalConstant * otherParticle.m_mass / distanceSquared;

            // Calculate the components of the gravitational acceleration
            alpha = std::asin(deltaY / std::sqrt(distanceSquared));
            float accelerationX = forceMagnitude * std::cos(alpha);
            float accelerationY = forceMagnitude * std::sin(alpha);

            // Add the components to the particle's acceleration
            particle.m_accel.x += accelerationX;
            particle.m_accel.y += accelerationY;
        }
    }
}

void Simulation::leapfrogUpdate(const float i_dt)
{
    // Update half-step velocities
    for (Particle& p : m_particleContainer)
    {
        p.m_vel_half_dt.x = p.m_vel.x + 0.5f * i_dt * p.m_accel.x;
        p.m_vel_half_dt.y = p.m_vel.y + 0.5f * i_dt * p.m_accel.y;
    }

    // Update positions
    for (Particle& p : m_particleContainer)
    {
        p.m_pos.x += i_dt * p.m_vel_half_dt.x;
        p.m_pos.y += i_dt * p.m_vel_half_dt.y;
    }

    // Calculate new accelerations
    for (Particle& p : m_particleContainer)
    {
        // Calculate new acceleration based on updated positions (can depend on external forces, etc.)
        calculateAcceleration(p);
    }

    // Update full-step velocities
    for (Particle& p : m_particleContainer)
    {
        p.m_vel.x = p.m_vel_half_dt.x + 0.5f * i_dt * p.m_accel.x;
        p.m_vel.y = p.m_vel_half_dt.y + 0.5f * i_dt * p.m_accel.y;
    }
}

void Simulation::moveParticles(const float i_dt)
{
    leapfrogUpdate(i_dt);
}

void Simulation::calcTotalEnergy()
{
    float kinEnergy_i{ 0.0F };
    float potEnergy_i{ 0.0F };
    float tmpVelScalarSqd{ 0.0F };
    float totalEnergy{ 0.0F };
    float dist{ 0.0F };
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        tmpVelScalarSqd = m_particleContainer[i].m_vel.x * m_particleContainer[i].m_vel.x + m_particleContainer[i].m_vel.y * m_particleContainer[i].m_vel.y;
        kinEnergy_i = 0.5 * m_particleContainer[i].m_mass * tmpVelScalarSqd;
        //sum of the pot energy with respect to all other particles
        for (size_t j = i+1; j < m_particleContainer.size(); j++)
        {
            dist = std::sqrt(std::pow((m_particleContainer[i].m_pos.x - m_particleContainer[j].m_pos.x), 2.0) + std::pow((m_particleContainer[i].m_pos.y - m_particleContainer[j].m_pos.y), 2.0));
            potEnergy_i += -m_particleContainer[i].m_mass * m_particleContainer[j].m_mass / dist;
        }
        totalEnergy += totalEnergy + potEnergy_i;
    }
    std::cout << totalEnergy << std::endl;
}