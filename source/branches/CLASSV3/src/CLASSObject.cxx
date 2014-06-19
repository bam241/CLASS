#include "CLASSObject.hxx"

#include "LogFile.hxx"

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

CLASSObject::CLASSObject(LogFile* log)
{
	SetLog(log);
}