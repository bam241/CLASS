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
}

CLASSFuelPlan::CLASSFuelPlan(CLASSLogger* log):CLASSObject(log)
{
}

void CLASSFuelPlan::AddFuel(cSecond time,  CLASSFuel fuel, double BurnUp)
{
	fLoadingPlan.insert(pair< cSecond, pair< CLASSFuel, double > >(time, pair< CLASSFuel, double > (fuel, BurnUp)) );
}

pair< CLASSFuel, double > CLASSFuelPlan::GetFuelAt(cSecond t)
{

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
			return FuelAtT;
		}
		else if (!AfterPrevious && !AfterNext)
		{
			WARNING << "The time asked is before the first laoding time..."<< endl;
			WARNING << "The first Fuel will be loaded... Check your FuelPLan!!!!!" << endl;
			return FuelAtT;
		}
		else
		{
			FuelAtT = (*it).second;
			AfterPrevious = true;
		}
	}

	return FuelAtT;
}
