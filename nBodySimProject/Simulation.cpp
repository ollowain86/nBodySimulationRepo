#define _USE_MATH_DEFINES //for M_PI
#include <cmath>
#include "Simulation.h"
#include "Particle.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <random>

//ctor
Simulation::Simulation(const size_t i_numOfParticles, const float i_scale, const float i_gravitationalConstant, const unsigned int i_edgeFreePixels) : m_numberOfParticles(i_numOfParticles), m_scale(i_scale), m_gravitationalConstant(i_gravitationalConstant), m_edgeFreePixels(i_edgeFreePixels)
{

}

// assigns pos, vel and mass to the particles
void Simulation::setUpSimulation(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr)
{
    /*
    Particle tmpParticle;
    for (size_t i = 0; i < m_numberOfParticles; i++)
    {
        if (i==0)
        {
            tmpParticle.m_mass = 10.0F;
            tmpParticle.m_pos.x = 1050;
            tmpParticle.m_pos.y = 500.0;
            tmpParticle.m_vel.x = 0.0;
            tmpParticle.m_vel.y = 2.0;
            tmpParticle.m_accel.x = 0.0F;
            tmpParticle.m_accel.y = 0.0F;
        }
        else
        {
            tmpParticle.m_mass = 1000.0F;
            tmpParticle.m_pos.x = 800.0;
            tmpParticle.m_pos.y = 500.0;
            tmpParticle.m_vel.x = 0.0;
            tmpParticle.m_vel.y = 0.0;
            tmpParticle.m_accel.x = 0.0F;
            tmpParticle.m_accel.y = 0.0F;
        }
        m_particleContainer.push_back(tmpParticle);
    }
    */
    
    //calculate plummer radius
    m_plummerRadiusSqd = (std::min(i_maxXlengthDistr, i_maxYlengthDistr) + std::abs(static_cast<float>(i_maxXlengthDistr - i_maxYlengthDistr)) / 2.0) / m_numberOfParticles;

    std::random_device rd;
    std::mt19937 gen(123456);
    double xDistrLengthMax = static_cast<double>(i_maxXlengthDistr);
    double yDistrLengthMax = static_cast<double>(i_maxYlengthDistr);
    float secondScale = m_scale/1.0;
    std::uniform_real_distribution<double> distr_x(static_cast<double>(m_edgeFreePixels) * secondScale, xDistrLengthMax - static_cast<double>(m_edgeFreePixels) * secondScale);
    std::uniform_real_distribution<double> distr_y(static_cast<double>(m_edgeFreePixels) * secondScale, yDistrLengthMax - static_cast<double>(m_edgeFreePixels) * secondScale);
    std::uniform_real_distribution<double> distrVel_x(-10.0, 10.0F);
    std::uniform_real_distribution<double> distrVel_y(-10.0, 10.0F);
    Particle tmpParticle;
    for (size_t i = 0; i < m_numberOfParticles; i++)
    {
        if (i == 0)
        {
            tmpParticle.m_mass = 10.0F;
        }
        else
        {
            tmpParticle.m_mass = 10.0F;
        }
        tmpParticle.m_pos.x = distr_x(gen);
        tmpParticle.m_pos.y = distr_y(gen);
        //tmpParticle.m_vel.x = distrVel_x(gen);
        //tmpParticle.m_vel.y = distrVel_y(gen);
        tmpParticle.m_vel.x = 0.0;
        tmpParticle.m_vel.y = 0.0;
        tmpParticle.m_accel.x = 0.0F;
        tmpParticle.m_accel.y = 0.0F;
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
    for (const Particle& otherParticle : m_particleContainer)
    {
        if (&particle != &otherParticle)  // Ensure not calculating for the same particle
        {
            // Calculate distance between the two particles
            m_deltaX = otherParticle.m_pos.x - particle.m_pos.x;
            m_deltaY = otherParticle.m_pos.y - particle.m_pos.y;
            m_distance = std::sqrt(m_deltaX * m_deltaX + m_deltaY * m_deltaY);

            // Calculate the gravitational force magnitude
            m_accelMagnitude = m_gravitationalConstant * otherParticle.m_mass / (m_distance*m_distance+m_plummerRadiusSqd);

            // Calculate the components of the gravitational acceleration
            m_accelerationX = m_accelMagnitude * (m_deltaX / m_distance);
            m_accelerationY = m_accelMagnitude * (m_deltaY / m_distance);

            // Add the components to the particle's acceleration
            particle.m_accel.x += m_accelerationX;
            particle.m_accel.y += m_accelerationY;
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
        p.m_vel_half_dt.x *= 1.0F;
        p.m_vel_half_dt.y *= 1.0F;
    }

    // Update positions
    for (Particle& p : m_particleContainer)
    {
        p.m_pos.x += i_dt * p.m_vel_half_dt.x;
        p.m_pos.y += i_dt * p.m_vel_half_dt.y;
    }

    // Calculate new accelerations O(N^2)
    for (Particle& p : m_particleContainer)
    {
        p.m_accel.x = 0.0F;
        p.m_accel.y = 0.0F;
        // Calculate new acceleration based on updated positions
        calculateAcceleration(p);
        p.m_accel.x *= 1.0F;
        p.m_accel.y *= 1.0F;
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
        totalEnergy += kinEnergy_i + potEnergy_i;
        potEnergy_i = 0.0F;
    }
    std::cout << totalEnergy << std::endl;
}

void Simulation::writeOutData()
{
    std::ofstream outFile("test.txt", std::ios::app);
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        outFile << m_particleContainer[i].m_pos.x << " " << m_particleContainer[i].m_pos.y << " " << std::sqrt((m_particleContainer[i].m_vel.x * m_particleContainer[i].m_vel.x) + (m_particleContainer[1].m_vel.y * m_particleContainer[1].m_vel.y)) << " " << std::sqrt((m_particleContainer[i].m_accel.x * m_particleContainer[i].m_accel.x) + (m_particleContainer[i].m_accel.y * m_particleContainer[i].m_accel.y));
        if (i == m_particleContainer.size() - 1)
        {
            outFile << "\n"; // Add a new line after every 4 particles or at the end of the loop
        }
        else
        {
            outFile << " "; // Add a space between each particle's data
        }
    }
    std::cout << std::endl;
    outFile.close();
}