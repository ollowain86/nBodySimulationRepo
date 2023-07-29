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
	//scale, gravitationalConstant, edgeFreePixels
	Simulation m_sim(2U, 1.0F, 1.0F, 0U);
	//maxXlengthDistr, maxYlengthDistr -> does not matter random values
	m_sim.setUpSelector(0U, 0U, 1U);
}