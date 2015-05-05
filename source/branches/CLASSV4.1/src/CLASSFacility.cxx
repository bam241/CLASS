#include "CLASSFacility.hxx"
#include "Scenario.hxx"

using namespace std;

	//________________________________________________________________________
	//
	//		CLASSFacility
	//
	//
	//
	//
	//________________________________________________________________________

ClassImp(CLASSFacility)




CLASSFacility::CLASSFacility(int type):CLASSObject()
{
	fParc = 0;
	
	fFacilityType = type;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = -1;
}

CLASSFacility::CLASSFacility(CLASSLogger* log, int type):CLASSObject(log)
{
	fParc = 0;

	fFacilityType = type;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = -1;
}


CLASSFacility::CLASSFacility(CLASSLogger* log, cSecond cycletime, int type):CLASSObject(log)
{
	fParc = 0;
	
	fFacilityType = type;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = cycletime;

}

CLASSFacility::CLASSFacility(CLASSLogger* log, cSecond creationtime, cSecond lifetime, int type):CLASSObject(log)
{
	fParc = 0;

	fFacilityType = type;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = -1;
	
	fCreationTime = creationtime;
	fLifeTime = lifetime;
}

CLASSFacility::CLASSFacility(CLASSLogger* log, cSecond creationtime, cSecond lifetime, cSecond cycletime, int type):CLASSObject(log)
{
	fParc = 0;

	fFacilityType = type;
	
	fInternalTime = fCreationTime;
	fInCycleTime = 0;
	fCycleTime = cycletime;
	
	fCreationTime = creationtime;
	fLifeTime = lifetime;
}



void CLASSFacility::ApplyZAIThreshold(int z)
{
	fInsideIV.ApplyZAIThreshold(z);
	fCumulativeIVIn.ApplyZAIThreshold(z);
	fCumulativeIVOut.ApplyZAIThreshold(z);

}