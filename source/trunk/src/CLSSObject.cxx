#include "CLSSObject.hxx"

#include "LogFile.hxx"

using namespace std;

	//________________________________________________________________________
	//
	//		CLSSObject
	//
	//
	//
	//
	//________________________________________________________________________

ClassImp(CLSSObject)



CLSSObject::CLSSObject()
{
}


CLSSObject::CLSSObject(LogFile* log)
{
	fLog = log;
}