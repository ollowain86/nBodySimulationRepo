#include <iostream>
#include "Visualizer.h"


int main()
{
    //number of particles
    const size_t numOfParticles{ 2U };
    const float scale = 1.0F;
    const float gravitationConstant = 1.0F;

    Visualizer visualizer(numOfParticles, scale, gravitationConstant);
    visualizer.render();

    return 0;
}