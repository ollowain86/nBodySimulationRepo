#include <iostream>
#include "Visualizer.h"
#include "TestinClass.h"


int main()
{
    //number of particles
    const size_t numOfParticles{ 500U };
    double scale = 1.0;
    const double gravitationConstant = 1.0;
    const unsigned int edgeFreePixels = 100U; //number of pixel to the edges of the monitor where no particle should be at start
    const unsigned int simOption = 0U; //0U random globular cluster, 1U two particle predefined test, 2U a milky way like structure

    //TestinClass tester;
    //tester.runTests();

    if (simOption == 1U)
    {
        scale = 1.0;
    }

    Visualizer visualizer(numOfParticles, scale, gravitationConstant, edgeFreePixels, simOption, 0.01);
    visualizer.render();

    return 0;
}