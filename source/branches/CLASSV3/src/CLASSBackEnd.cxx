#include "CLASSBackEnd.hxx"

#include "DecayDataBank.hxx"
#include "CLASS.hxx"
#include "LogFile.hxx"


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

CLASSBackEnd::CLASSBackEnd():CLASSFacility()
{
	fDecayDataBase = 0;
}

CLASSBackEnd::CLASSBackEnd(LogFile* log):CLASSFacility(log)
{
	fDecayDataBase = 0;
}


CLASSBackEnd::CLASSBackEnd(LogFile* log, cSecond cycletime):CLASSFacility(log, cycletime)
{
	fDecayDataBase = 0;
}

//________________________________________________________________________
void CLASSBackEnd::ClearIVArray()
{

	IsotopicVector EmptyIV;
	fInsideIV = EmptyIV;
	fIVArray.clear();
}

//________________________________________________________________________
void CLASSBackEnd::AddIV(IsotopicVector isotopicvector)
{

	AddCumulativeIVIn(isotopicvector);

	fIVArray.push_back(isotopicvector);

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
