#include <iostream>
#include "Visualizer.h"
#include "TestinClass.h"


int main()
{
    //number of particles
    const size_t numOfParticles{ 2U };
    const float scale = 10.0F;
    const float gravitationConstant = 1.0F;
    const unsigned int edgeFreePixels = 100U; //number of pixel to the edges of the monitor where no particle should be at start

    //TestinClass tester;
    //tester.runTests();

    Visualizer visualizer(numOfParticles, scale, gravitationConstant, edgeFreePixels);
    visualizer.render();

    return 0;
}