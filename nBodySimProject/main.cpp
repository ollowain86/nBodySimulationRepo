#include <iostream>
#include "Visualizer.h"


int main()
{
    //number of particles
    const size_t numOfParticles{ 2U };
    const float scale = 1.0F;
    const float gravitationConstant = 1.0F;
    const unsigned int edgeFreePixels = 50U; //number of pixel to the edges of the monitor where no particle should be at start

    Visualizer visualizer(numOfParticles, scale, gravitationConstant, edgeFreePixels);
    visualizer.render();

    return 0;
}