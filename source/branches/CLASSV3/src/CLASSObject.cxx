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
	fLog = new CLASSLogger("CLASSObject.log");
}

CLASSObject::CLASSObject(CLASSLogger* log)
{
	fLog = log;
}