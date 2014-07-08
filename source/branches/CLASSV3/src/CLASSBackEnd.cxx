#include "CLASSBackEnd.hxx"

#include "DecayDataBank.hxx"
#include "Scenario.hxx"
#include "CLASSLogger.hxx"


#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

//________________________________________________________________________
//
//		CLASSBackEnd
//
//
//
//
//________________________________________________________________________
ClassImp(CLASSBackEnd)

CLASSBackEnd::CLASSBackEnd(int type):CLASSFacility(type)
{
	fDecayDataBase = 0;
}

CLASSBackEnd::CLASSBackEnd(CLASSLogger* log, int type):CLASSFacility(log, type)
{
	fDecayDataBase = 0;
}


CLASSBackEnd::CLASSBackEnd(CLASSLogger* log, cSecond cycletime, int type):CLASSFacility(log, cycletime, type)
{
	fDecayDataBase = 0;
}

//________________________________________________________________________
void CLASSBackEnd::ClearIVArray()
{

	IsotopicVector EmptyIV;
	fInsideIV = EmptyIV;
	fIVArray.clear();
	fIVArrayArrivalTime.clear();
}

//________________________________________________________________________
void CLASSBackEnd::AddIV(IsotopicVector isotopicvector)
{

	AddCumulativeIVIn(isotopicvector);

	fIVArray.push_back(isotopicvector);
	fIVArrayArrivalTime.push_back(fInternalTime);


}

//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector CLASSBackEnd::GetDecay(IsotopicVector isotopicvector, cSecond t)
{

	IsotopicVector IV;

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->Evolution(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}

	return IV;

}
