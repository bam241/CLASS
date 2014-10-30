#include "CLASSFuelPlan.hxx"

#include "CLASSLogger.hxx"

using namespace std;

//________________________________________________________________________
//
//		CLASSFuelPlan
//
//
//
//
//________________________________________________________________________

CLASSFuelPlan::CLASSFuelPlan():CLASSObject()
{
	DBGL
}

CLASSFuelPlan::CLASSFuelPlan(CLASSLogger* log):CLASSObject(log)
{
	DBGL
	DBGL
}

void CLASSFuelPlan::AddFuel(cSecond time,  CLASSFuel fuel, double BurnUp)
{
	DBGL
	fLoadingPlan.insert(pair< cSecond, pair< CLASSFuel, double > >(time, pair< CLASSFuel, double > (fuel, BurnUp)) );
	DBGL
}

pair< CLASSFuel, double > CLASSFuelPlan::GetFuelAt(cSecond t)
{
	DBGL

	pair< CLASSFuel, double > FuelAtT = fLoadingPlan.begin()->second;

	map< cSecond, pair< CLASSFuel, double > >::iterator it;

	cSecond ThisFuelTime;

	for (it = fLoadingPlan.begin(); it != fLoadingPlan.end(); it++ )
	{
		if( it == fLoadingPlan.begin())
		{
			FuelAtT = (*it).second;
		}
		else
		{
			ThisFuelTime = (*it).first;

			if( t < ThisFuelTime )
			{
				DBGL
				return FuelAtT;
			}
			else
			{
				FuelAtT = (*it).second;
			}

		}
	}

	DBGL
	return FuelAtT;
}
