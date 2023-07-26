#include <iostream>
#include "Visualizer.h"


int main()
{
    //number of particles
    size_t numOfParticles{ 2000U };
    float scale = 1000.0F;

    Visualizer visualizer(numOfParticles, scale);
    visualizer.render();

    return 0;
}