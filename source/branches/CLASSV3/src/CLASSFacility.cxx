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


CLASSFacility::CLASSFacility():CLASSObject()
{
	fParc = 0;
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = -1;
}

CLASSFacility::CLASSFacility(LogFile* log):CLASSObject(log)
{
	fParc = 0;
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = -1;
}


CLASSFacility::CLASSFacility(LogFile* log, cSecond cycletime):CLASSObject(log)
{
	fParc = 0;
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = cycletime;
}

CLASSFacility::CLASSFacility(LogFile* log, cSecond creationtime, cSecond lifetime):CLASSObject(log)
{
	fParc = 0;
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = -1;
	fCreationTime = creationtime;
	fLifeTime = lifetime;
}

CLASSFacility::CLASSFacility(LogFile* log, cSecond creationtime, cSecond lifetime, cSecond cycletime):CLASSObject(log)
{
	fParc = 0;
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = cycletime;
	fCreationTime = creationtime;
	fInternalTime = fCreationTime;
	fLifeTime = lifetime;
}