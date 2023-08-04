#include "TestinClass.h"
#include "Simulation.h"

TestinClass::TestinClass()
{

}

void TestinClass::runTests()
{
	test1();
}

void TestinClass::test1()
{
	//´num of particles, scale, gravitationalConstant, edgeFreePixels, option
	Simulation m_sim(2U, 1.0, 1.0, 0U, 1U);
	//maxXlengthDistr, maxYlengthDistr -> does not matter random values
	m_sim.setUpSelector(0U, 0U);
}