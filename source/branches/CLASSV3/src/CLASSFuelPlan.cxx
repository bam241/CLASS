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
	bool AfterPrevious = false;
	bool AfterNext = true;
	for (it = fLoadingPlan.begin(); it != fLoadingPlan.end(); it++ )
	{
		if (t < (*it).first )
			AfterNext = false;
		else
			AfterNext = true;

		if (AfterPrevious && !AfterNext)
		{
			DBGL
			return FuelAtT;
		}
		else if (!AfterPrevious && !AfterNext)
		{
			WARNING << "The time asked is before the first laoding time..."<< endl;
			WARNING << "The first Fuel will be loaded... Check your FuelPLan!!!!!" << endl;
			DBGL
			return FuelAtT;
		}
		else
		{
			FuelAtT = (*it).second;
			AfterPrevious = true;
		}
	}

	DBGL
	return FuelAtT;
}
