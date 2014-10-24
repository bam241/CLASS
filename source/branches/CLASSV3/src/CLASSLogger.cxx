
#include "CLASSLogger.hxx"

#include <iostream>
#include <ostream>
#include <algorithm>
using namespace std;

//________________________________________________________________________
//
//		CLASSLogger
//
//
//
//
//________________________________________________________________________
CLASSLogger::CLASSLogger()
{
	string CLASSLoggerName = "CLASS_OUTPUT.log";
	int VerboseLvl = 0;
	int OutputLvl = 1;

	
	fCLASSLoggerName = CLASSLoggerName;
	fOutPutFile.open(CLASSLoggerName.c_str());
	fVerboseLVL = VerboseLvl;
	if(!fOutPutFile)
	{
		cout << "Could not open the CLASSLogger: " << CLASSLoggerName << " !" << endl;
		exit(-1);
	}

	fError = 0;
	fWarning = 0;
	fDebug = 0;
	fInfo = 0;

	if (VerboseLvl <= OutputLvl)
		fMaxOutPutLVL = OutputLvl;
	else
		fMaxOutPutLVL = VerboseLvl;

	if(VerboseLvl >= 0)
	{
		if (!fError)
			fError = new LogType(std::cout);
		else
			fError->SetSecondOutput(std::cout);
	}
	if(VerboseLvl >= 1)
	{
		if (!fWarning)
			fWarning = new LogType(std::cout);
		else
			fWarning->SetSecondOutput(std::cout);
	}
	if(VerboseLvl >= 2)
	{
		if (!fInfo)
			fInfo = new LogType(std::cout);
		else
			fInfo->SetSecondOutput(std::cout);
	}
	if(VerboseLvl >= 3)
	{
		if (!fDebug)
			fDebug = new LogType(std::cout);
		else
			fDebug->SetSecondOutput(std::cout);
	}
	if(OutputLvl >= 0)
	{
		if (!fError)
			fError = new LogType(fOutPutFile);
		else
			fError->SetSecondOutput(fOutPutFile);
	}
	if(OutputLvl >= 1)
	{
		if (!fWarning)
			fWarning = new LogType(fOutPutFile);
		else
			fWarning->SetSecondOutput(fOutPutFile);
	}
	if(OutputLvl >= 2)
	{
		if (!fInfo)
			fInfo = new LogType(fOutPutFile);
		else
			fInfo->SetSecondOutput(fOutPutFile);
	}
	if(OutputLvl >= 3)
	{
		if (!fDebug)
			fDebug = new LogType(fOutPutFile);
		else
			fDebug->SetSecondOutput(fOutPutFile);
	}
	
	
	
}



CLASSLogger::CLASSLogger(string CLASSLoggerName, int VerboseLvl, int OutputLvl )
{
	fCLASSLoggerName = CLASSLoggerName;
	fOutPutFile.open(CLASSLoggerName.c_str());
	fVerboseLVL = VerboseLvl; 
	if(!fOutPutFile)
	{
		cout << "Could not open the CLASSLogger: " << CLASSLoggerName << " !" << endl;
		exit(-1);
	}

	fError = 0;
	fWarning = 0;
	fDebug = 0;
	fInfo = 0;

	if (VerboseLvl <= OutputLvl)
		fMaxOutPutLVL = OutputLvl;
	else
		fMaxOutPutLVL = VerboseLvl;

	if(VerboseLvl >= 0)
	{
		if (!fError)
			fError = new LogType(std::cout);
		else
			fError->SetSecondOutput(std::cout);
	}

	if(VerboseLvl >= 1)
	{
		if (!fWarning)
			fWarning = new LogType(std::cout);
		else
			fWarning->SetSecondOutput(std::cout);
	}
	if(VerboseLvl >= 2)
	{
		if (!fInfo)
			fInfo = new LogType(std::cout);
		else
			fInfo->SetSecondOutput(std::cout);
	}
	if(VerboseLvl >= 3)
	{
	if (!fDebug)
			fDebug = new LogType(std::cout);
		else
			fDebug->SetSecondOutput(std::cout);
	}


	if(OutputLvl >= 0)
	{
		if (!fError)
			fError = new LogType(fOutPutFile);
		else
			fError->SetSecondOutput(fOutPutFile);
	}
	if(OutputLvl >= 1)
	{
		if (!fWarning)
			fWarning = new LogType(fOutPutFile);
		else
			fWarning->SetSecondOutput(fOutPutFile);
	}
	if(OutputLvl >= 2)
	{
		if (!fInfo)
			fInfo = new LogType(fOutPutFile);
		else
			fInfo->SetSecondOutput(fOutPutFile);
	}
	if(OutputLvl >= 3)
	{
		if (!fDebug)
			fDebug = new LogType(fOutPutFile);
		else
			fDebug->SetSecondOutput(fOutPutFile);
	}



}

//________________________________________________________________________
CLASSLogger::~CLASSLogger()
{

	fOutPutFile.close();

}




