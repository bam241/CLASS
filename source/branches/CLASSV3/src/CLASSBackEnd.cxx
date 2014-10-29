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
void CLASSBackEnd::UpdateInsideIV()
{
	DBGL
	fInsideIV = IsotopicVector();
	for(int i = 0; i < (int)fIVArray.size(); i++)
		fInsideIV += fIVArray[i];
	DBGL
}


//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector CLASSBackEnd::GetDecay(IsotopicVector isotopicvector, cSecond t)
{
	DBGL

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

	DBGL
	return IV;
}


map<cSecond,int> CLASSBackEnd::GetTheBackEndTimePath()
{
	DBGL
	map<cSecond, int> TheBackEndTimePath;

	if(!fIsStorageType)
	{
		int FacilityType = GetFacilityType();
		cSecond step = GetCycleTime();

		pair< map<cSecond, int>::iterator, bool > IResult  = TheBackEndTimePath.insert(pair<cSecond,int> (step, FacilityType));
		if( !IResult.second ) IResult.first->second |= FacilityType;

		map<cSecond, int> TheBackEndTimePath_tmp = GetOutBackEndFacility()->GetTheBackEndTimePath();

		map<cSecond, int>::iterator it;
		for (it = TheBackEndTimePath_tmp.begin(); it != TheBackEndTimePath_tmp.end(); it++)
		{
			pair< map<cSecond, int>::iterator, bool > IResult;

			IResult = TheBackEndTimePath.insert( pair<cSecond ,int>(step + (*it).first, (*it).second) );
			if( !IResult.second )
				IResult.first->second |= (*it).second;

		}
	}

	DBGL
	return TheBackEndTimePath;
}