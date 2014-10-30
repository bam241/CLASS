#include "SeparationPlant.hxx"

#include "IsotopicVector.hxx"
#include "Storage.hxx"
#include "Scenario.hxx"
#include "CLASSLogger.hxx"

#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

//________________________________________________________________________
//
//		SeparationPlant
//
//
//
//
//________________________________________________________________________
ClassImp(SeparationPlant)


SeparationPlant::SeparationPlant():CLASSBackEnd(-2)
{
	fOutBackEndFacility = 0;
	SetName("C_SeparationPlant.");
	SetIsStorageType();
}

//________________________________________________________________________
SeparationPlant::SeparationPlant(CLASSLogger* log):CLASSBackEnd(log, -2)
{


	fCycleTime = 0;
	fIsStarted = false;
	fPutToWaste = true;
	fCoolingLastIndex = 0;

	fOutBackEndFacility = 0;
	SetName("C_SeparationPlant.");

	SetIsStorageType();
	
	INFO	<< " A new SeparationPlant has been define :" << endl;
	INFO	<< "\t The Separation Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	WARNING	<< " All Separated Fuel go directly to WASTE after cooling !! " << endl;


}


//________________________________________________________________________
SeparationPlant::~SeparationPlant()
{


}

//________________________________________________________________________

//________________________________________________________________________
void SeparationPlant::SetBackEndDestination(CLASSBackEnd* storagedestination, IsotopicVector isotopicvector, cSecond destinationstartingtime)
{
	DBGL

	fDestinationStorageStartingTime.push_back(destinationstartingtime);
	fDestinationStorage.push_back(storagedestination);
	fDestinationStorageIV.push_back(isotopicvector);

	if(fDestinationStorage.size()>=2)
	{
		//checker que pas 2 fois le même ZAI
		//Dans VectorIsotopic, BaM a fait des trucs cools
		//-> IV->GetThisComposition(IV2)
	}
	if (fDestinationStorage.size() != fDestinationStorageIV.size())
	ERROR	<< " fDestinationStorage.size() != fDestinationStorageIV.size() !! " << endl;


/*



les pertes non gérées -> WASTE par défaut


	

*/
	DBGL

}


void SeparationPlant::AddIV(IsotopicVector IV)
{ 
DBGL
	for(int fds=0; fds<(int)fDestinationStorage.size(); fds++)
	{
		cSecond CurrentTime = GetParc()->GetAbsoluteTime();

		INFO << "Separation..." <<endl;
		INFO << "Current Time : " << CurrentTime <<endl;
		INFO << "IV Separation Time : " << fDestinationStorageStartingTime[fds] <<endl;

		if(CurrentTime >= fDestinationStorageStartingTime[fds])
		{
			IsotopicVector IVtmp;
			IVtmp = IV*fDestinationStorageIV[fds];
			fDestinationStorage[fds]->AddIV(IVtmp);
			IV -= IVtmp;
		}

	}

	GetParc()->AddWaste(IV);

	DBGL
}



map<cSecond,int> SeparationPlant::GetTheBackEndTimePath()
{
	DBGL

	map<cSecond, int> TheBackEndTimePath;
	for( int i = 0; i < (int)fDestinationStorage.size(); i++ )
	{
		map<cSecond, int> TheBackEndTimePath_tmp = fDestinationStorage[i]->GetTheBackEndTimePath();
		map<cSecond, int>::iterator it;
		for (it = TheBackEndTimePath_tmp.begin(); it != TheBackEndTimePath_tmp.end(); it++)
		{
			pair< map<cSecond, int>::iterator, bool > IResult;
			IResult = TheBackEndTimePath.insert( pair<cSecond ,int>((*it).first, (*it).second) );
			if( !IResult.second )
				IResult.first->second |= (*it).second;

		}
	}

	DBGL
	return TheBackEndTimePath;
}





//________________________________________________________________________
//	Time Action with the reste of the Universe : 
//		Out Storage
//		Evolution : 
//			Cooling
//________________________________________________________________________



