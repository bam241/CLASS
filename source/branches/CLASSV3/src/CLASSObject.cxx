#include "CLASSObject.hxx"

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
}

CLASSObject::CLASSObject(CLASSLogger* log)
{
	fLog = log;
}