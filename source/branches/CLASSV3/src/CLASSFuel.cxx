#include "CLASSFuel.hxx"

#include "CLASSLogger.hxx"

using namespace std;

//________________________________________________________________________
//
//		CLASSFuel
//
//
//
//
//________________________________________________________________________

CLASSFuel::CLASSFuel(EvolutionData* evo)
{
	fEvolutionData = evo;
	fPhysicsModels = 0;
}

CLASSFuel::CLASSFuel(PhysicsModels* evo)
{
	fEvolutionData = 0;
	fPhysicsModels = evo;
}