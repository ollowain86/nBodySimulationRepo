#define _USE_MATH_DEFINES //for M_PI
#include <cmath>
#include "Simulation.h"
#include "Particle.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include "helperTypesAndFunctions.h"

//ctor
Simulation::Simulation(const size_t i_numOfParticles, const double i_scale, const double i_gravitationalConstant, const unsigned int i_edgeFreePixels) : m_numberOfParticles(i_numOfParticles), m_scale(i_scale), m_gravitationalConstant(i_gravitationalConstant), m_edgeFreePixels(i_edgeFreePixels)
{

}

void Simulation::setUpSelector(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr, const unsigned int i_option)
{
    // globular cluster circular shape - virialized
    if (i_option == 0U)
    {
        setUpSimulation(i_maxXlengthDistr, i_maxYlengthDistr);
    }
    // hardcoded 2 particle system
    else if(i_option == 1U)
    {
        setUpTwoParticle();
    }
    // hardcoded 3 particle system
    else if (i_option == 2U)
    {

    }
    else if (i_option == 3U)
    {

    }
    else
    {
        std::cout << "In setUpSelector, not-defined option is selected." << std::endl;
    }
}

bool Simulation::doubleEqual(const double a, const double b)
{
    bool isEqual = false;
    double epsilon = 10e-5;

    if (std::abs(a-b)<epsilon)
    {
        isEqual = true;
    }
    return isEqual;
}

void Simulation::setUpTwoParticle()
{
    // Circular orbit: v = sqrt(G*M/radius)
    Particle tmpParticle;
    for (size_t i = 0; i < 2U; i++)
    {
        if (i == 0)
        {
            tmpParticle.m_mass = 1000.0;
            tmpParticle.m_pos.x = 500.0;
            tmpParticle.m_pos.y = 500.0;
            tmpParticle.m_vel.x = 0.0;
            tmpParticle.m_vel.y = 0.0;
            tmpParticle.m_accel.x = 0.0;
            tmpParticle.m_accel.y = 0.0;
        }
        else
        {
            tmpParticle.m_mass = 10.0;
            tmpParticle.m_pos.x = 510.0;
            tmpParticle.m_pos.y = 500.0;
            tmpParticle.m_vel.x = 0.0;
            tmpParticle.m_vel.y = 10.0;
            tmpParticle.m_accel.x = 0.0;
            tmpParticle.m_accel.y = 0.0;
        }
        m_particleContainer.push_back(tmpParticle);
    }
    // should not matter, since for a radius > 2*plummer it is normal newton
    m_plummerRadius = 10.0/2U;

    // calc accel initially
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        calculateAcceleration(m_particleContainer[i]);
    }

    double expected_ax{0.0};
    double expected_ay{ 0.0 };
    double expected_a{ 0.0 };
    double distance{ 0.0 };
    double deltaX{ 0.0 };
    double deltaY{ 0.0 };
    double accelMagnitude{ 0.0 };
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        for (size_t j = 0; j < m_particleContainer.size(); j++)
        {
            if (i != j)
            {
                deltaX = m_particleContainer[j].m_pos.x - m_particleContainer[i].m_pos.x;
                deltaY = m_particleContainer[j].m_pos.y - m_particleContainer[i].m_pos.y;
                distance = std::sqrt(deltaX * deltaX + deltaY * deltaY);
                accelMagnitude += m_particleContainer[j].m_mass/(distance* distance);
                if (i == 0U)
                {
                    //calculated manually
                    expected_a = 1/10.0;
                    expected_ax = expected_a * 10.0/10.0;
                    expected_ay = expected_a * 0.0 / 10.0;
                }
                else
                {
                    //calculated manually
                    expected_a = 10.0;
                    expected_ax = -10.0;
                    expected_ay = 0.0;
                }
                if (!doubleEqual(accelMagnitude, expected_a))
                {
                    std::cout << "Two Particle Sys: accelMagnitude wrong!" << std::endl;
                }
                if (!doubleEqual(accelMagnitude*deltaX/distance, expected_ax) || !doubleEqual(accelMagnitude * deltaY / distance, expected_ay))
                {
                    std::cout << "Two Particle Sys: ax or ay is wrong!" << std::endl;
                }
            }
            accelMagnitude = 0.0;
        }
    }
}

// particles positions are set randomly into a circular shape with their radius (no vel, no accel assigned here), but m_particleContainer gets its initial size
void Simulation::setUpCircularShape(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr, const double i_maxRadius)
{
    std::random_device rd;
    std::mt19937 gen(123456);

    // set particles around center of visualization -> calculate center
    double originX = static_cast<double>(i_maxXlengthDistr) / 2.0;
    double originY = static_cast<double>(i_maxYlengthDistr) / 2.0;
    double tmpRadius{ 0.0 };
    double tmpAngle{ 0.0 };

    //distribute particle with random radius (0 and maxRadius) and angle 0, 2*pi
    std::uniform_real_distribution<double> distr_radius(0.0, i_maxRadius);
    std::uniform_real_distribution<double> distr_angle(0.0, 2.0 * M_PI);

    Particle tmpParticle;
    for (size_t i = 0; i < m_numberOfParticles; i++)
    {
        tmpParticle.m_mass = 10.0;

        tmpRadius = distr_radius(gen);
        tmpAngle = distr_angle(gen);

        tmpParticle.m_radius = tmpRadius;
        tmpParticle.m_pos.x = originX + (tmpRadius * cos(tmpAngle));
        tmpParticle.m_pos.y = originY + (-1.0 * tmpRadius * sin(tmpAngle));

        tmpParticle.m_pos.x *= m_scale;
        tmpParticle.m_pos.y *= m_scale;

        m_particleContainer.push_back(tmpParticle);
    }
}

//calculates M_i (mass inside r_i)
double Simulation::massWithinRadiusCalculator(const Particle& i_particle)
{
    double xPos{ 0.0 };
    double yPos{ 0.0 };
    double massWithinR_i{ 0.0 };

    for (const Particle& otherparticle : m_particleContainer)
    {
        if (&otherparticle != &i_particle)
        {
            if (otherparticle.m_radius <= i_particle.m_radius)
            {
                massWithinR_i += otherparticle.m_mass;
            }
        }
    }

    return massWithinR_i;
}

// calcDistance between two particles
double Simulation::calcDistance(const Particle& particleA, const Particle& particleB)
{
    return std::sqrt(std::pow((particleA.m_pos.x-particleB.m_pos.x),2.0) + std::pow((particleA.m_pos.y - particleB.m_pos.y), 2.0));
}

double Simulation::calcTotalPotentialEnergy()
{
    double U{ 0.0 };

    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        for (size_t j = i+1; j < m_particleContainer.size(); j++)
        {
            U += -m_gravitationalConstant*m_particleContainer[i].m_mass * m_particleContainer[j].m_mass / calcDistance(m_particleContainer[i], m_particleContainer[j]);
        }
    }
    return U;
}

// calculate length of vec2f
double Simulation::calcLength(const helpers::vec2d_double i_2dVec)
{
    return std::sqrt(i_2dVec.x* i_2dVec.x + i_2dVec.y* i_2dVec.y);
}

double Simulation::calcTotalKineticEnergy()
{
    double T{ 0.0 };

    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        T += m_particleContainer[i].m_mass * m_particleContainer[i].m_velScalar * m_particleContainer[i].m_velScalar;
    }
    T *= 0.5F;
    return T;
}

void Simulation::calcOrbitalSpeed(Particle& particle, const double i_massWithinR_i)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    if (doubleEqual(particle.m_radius, 0.0))
    {
        particle.m_radius += 0.001F;
    }
    double velScalar = std::sqrt(m_gravitationalConstant*i_massWithinR_i/ particle.m_radius);
    std::uniform_real_distribution<double> distr_vel(-velScalar, velScalar);

    particle.m_vel.x = distr_vel(gen);
    
    particle.m_vel.y = std::sqrt(velScalar*velScalar - particle.m_vel.x* particle.m_vel.x);
    
    particle.m_velScalar = velScalar;
}

// assigns pos, vel and mass to the particles
void Simulation::setUpSimulation(const unsigned int i_maxXlengthDistr, const unsigned int i_maxYlengthDistr)
{   
    double maxRadius = (static_cast<double>(std::min(i_maxXlengthDistr, i_maxYlengthDistr)) - 2.0*static_cast<double>(m_edgeFreePixels)) / 2.0;
    // Calculate plummer radius = radius of simulation / number of particles
    m_plummerRadius = maxRadius / m_numberOfParticles;

    // particles positions are set randomly into a circular shape with their radius (no vel, no accel assigned here), but m_particleContainer gets its initial size
    setUpCircularShape(i_maxXlengthDistr, i_maxYlengthDistr, maxRadius);

    double massWithinR_i{ 0.0 };
    //first calculate the mass for particle with mass m_i inside its radius r_i
    //second calculate the orbital velocity with the mass
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        //calculates M_i (mass inside r_i)
        massWithinR_i = massWithinRadiusCalculator(m_particleContainer[i]);
        //calc initial orbital speed
        if (massWithinR_i == 9990.0)
        {
            calcOrbitalSpeed(m_particleContainer[i], massWithinR_i);
        }
        calcOrbitalSpeed(m_particleContainer[i], massWithinR_i);
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

double Simulation::calcTotalEnergy()
{
    double kinEnergy_i{ 0.0 };
    double potEnergy_i{ 0.0 };
    double tmpVelScalarSqd{ 0.0 };
    double totalEnergy{ 0.0 };
    double dist{ 0.0 };
    for (size_t i = 0; i < m_particleContainer.size(); i++)
    {
        tmpVelScalarSqd = m_particleContainer[i].m_vel.x * m_particleContainer[i].m_vel.x + m_particleContainer[i].m_vel.y * m_particleContainer[i].m_vel.y;
        kinEnergy_i = 0.5 * m_particleContainer[i].m_mass * tmpVelScalarSqd;
        //sum of the pot energy with respect to all other particles
        for (size_t j = i + 1; j < m_particleContainer.size(); j++)
        {
            dist = std::sqrt(std::pow((m_particleContainer[i].m_pos.x - m_particleContainer[j].m_pos.x), 2.0) + std::pow((m_particleContainer[i].m_pos.y - m_particleContainer[j].m_pos.y), 2.0));
            potEnergy_i += -m_particleContainer[i].m_mass * m_particleContainer[j].m_mass / dist;
        }
        totalEnergy += kinEnergy_i + potEnergy_i;
        potEnergy_i = 0.0;
    }
    return totalEnergy;
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
            m_distanceSqrd = m_distance * m_distance;

            // Calculate the gravitational force magnitude
            if (m_distance < 2.0*m_plummerRadius)
            {
                m_accelMagnitude = m_gravitationalConstant * otherParticle.m_mass * (64.0 * m_distance * std::pow(m_plummerRadius, 3.0)) / std::pow((m_distanceSqrd + 4.0 * m_plummerRadius * m_plummerRadius), 3.0);
            }
            else
            {
                m_accelMagnitude = m_gravitationalConstant * otherParticle.m_mass / m_distanceSqrd;
            }

            // Calculate the components of the gravitational acceleration
            m_accelerationX = m_accelMagnitude * (m_deltaX / m_distance);
            m_accelerationY = m_accelMagnitude * (m_deltaY / m_distance);

            // Add the components to the particle's acceleration
            particle.m_accel.x += m_accelerationX;
            particle.m_accel.y += m_accelerationY;
        }
    }
}

void Simulation::leapfrogUpdate(const double i_dt)
{
    // calc r_(i+1/2)
    for (Particle& p : m_particleContainer)
    {
        p.m_pos.x = p.m_pos.x + p.m_vel.x * 0.5 * i_dt;
        p.m_pos.y = p.m_pos.y + p.m_vel.y * 0.5 * i_dt;
    }
    // calc a_(i+1/2)
    for (Particle& p : m_particleContainer)
    {
        // Calculate new acceleration based on updated positions
        p.m_accel.x = 0.0;
        p.m_accel.y = 0.0;
        calculateAcceleration(p);
    }
    // calc v_(i+1) with a_(i+1/2)
    for (Particle& p : m_particleContainer)
    {
        p.m_vel.x = p.m_vel.x + p.m_accel.x * i_dt;
        p.m_vel.y = p.m_vel.y + p.m_accel.y * i_dt;
        p.m_velScalar = calcLength(p.m_vel);
    }

    // calc r_(i+1) with v_(i+1)
    for (Particle& p : m_particleContainer)
    {
        p.m_pos.x = p.m_pos.x + p.m_vel.x * 0.5 * i_dt;
        p.m_pos.y = p.m_pos.y + p.m_vel.y * 0.5 * i_dt;
    }
}

void Simulation::moveParticles(const double i_dt)
{
    leapfrogUpdate(i_dt);
}

void Simulation::writeOutData()
{
    double T{ 0.0 };
    double U{ 0.0 };
    double E_tot{ 0.0 };
    std::ofstream outFile("test.txt", std::ios::app);
    T = calcTotalKineticEnergy();
    U = calcTotalPotentialEnergy();
    E_tot = calcTotalEnergy();
    outFile << T << " " << U << " " << T + U << " " << E_tot << std::endl;
    outFile.close();
}