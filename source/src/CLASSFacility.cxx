#include "CLASSFacility.hxx"
#include "Scenario.hxx"

	//________________________________________________________________________
	//
	//		CLASSFacility
	//
	//
	//
	//
	//________________________________________________________________________

ClassImp(CLASSFacility)

//____________________________________________________________________________
CLASSFacility::CLASSFacility ( int type ) :
	CLASSObject()    , fInternalTime(0) , fInCycleTime(0) , fCycleTime(-1) , fFacilityType(type) , fParc(nullptr)
{ ; }
//____________________________________________________________________________
CLASSFacility::CLASSFacility ( CLASSLogger* log , int type ) :
	CLASSObject(log) , fInternalTime(0) , fInCycleTime(0) , fCycleTime(-1) , fFacilityType(type) , fParc(nullptr)
{ ; }

//____________________________________________________________________________
CLASSFacility::CLASSFacility ( CLASSLogger* log , cSecond cycletime, int type) :
	CLASSObject(log) , fInternalTime(0) , fInCycleTime(0) , fCycleTime(cycletime) , fFacilityType(type) , fParc(nullptr)
{ ; }
//____________________________________________________________________________
CLASSFacility::CLASSFacility ( CLASSLogger* log , cSecond creationtime , cSecond lifetime , int type ) :
	CLASSObject(log) , fInternalTime(creationtime) , fInCycleTime(0) , fCycleTime(-1) , fFacilityType(type) , fParc(nullptr) , fCreationTime(creationtime) , fLifeTime(lifetime)
{ ; }
//____________________________________________________________________________
CLASSFacility::CLASSFacility(CLASSLogger* log, cSecond creationtime, cSecond lifetime, cSecond cycletime, int type):CLASSObject(log) , fInternalTime(creationtime) , fInCycleTime(0) , fCycleTime(cycletime) , fFacilityType(type) , fParc(nullptr) , fCreationTime(creationtime) , fLifeTime(lifetime)
{ ; }

//____________________________________________________________________________
void CLASSFacility::ApplyZAIThreshold(int z)
{
	fInsideIV.ApplyZAIThreshold(z);
	fCumulativeIVIn.ApplyZAIThreshold(z);
	fCumulativeIVOut.ApplyZAIThreshold(z);
}
