#include <iostream>
#include "Visualizer.h"


int main()
{
    //number of particles
    size_t numOfParticles{ 3U };

    Visualizer visualizer(numOfParticles);
    visualizer.render();

    return 0;
}