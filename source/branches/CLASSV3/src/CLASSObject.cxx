#include "CLASSObject.hxx"

#include "CLASSLogger.hxx"

using namespace std;

	//________________________________________________________________________
	//
	//		CLASSObject
	//
	//
	//
	//
	//________________________________________________________________________

ClassImp(CLASSObject)

CLASSObject::CLASSObject()
{
	fLog = 0;
	fIsLog = false;
}

CLASSObject::CLASSObject(CLASSLogger* log)
{
	SetLog(log);
}