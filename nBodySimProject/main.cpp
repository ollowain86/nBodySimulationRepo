#include <iostream>
#include "SFML/Graphics.hpp"
#include "Simulation.h"
#include <vector>
#include "Particle.h"

int main()
{
    //number of particles
    size_t numOfParticles{ 2U };
    Simulation simulation(numOfParticles);


    std::vector<Particle> particleContainer;

    Particle tmpParticle;
    for (size_t i = 0; i < numOfParticles; i++)
    {
        tmpParticle.setMass(10.0F);
        tmpParticle.setPos({i*10.0F, i*10.0F});
        tmpParticle.setPos({ i*3.0F, i*3.0F });

        particleContainer.push_back(tmpParticle);
    }

    return 0;
}