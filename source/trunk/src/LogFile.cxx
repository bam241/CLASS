
#include "LogFile.hxx"

#include <iostream>
#include <algorithm>
using namespace std;

//________________________________________________________________________
//
//		LogFile
//
//
//
//
//________________________________________________________________________

LogFile::LogFile(string LogFileName )
{
	fLogFileName = LogFileName;

	fLog.open(LogFileName.c_str());
	if(!fLog)
	{
		cout << "Could not open the LogFile: " << LogFileName << " !" << endl;
		exit(-1);
	}
	else
		cout << "LogFile: " << LogFileName << " opened." << endl;

}

//________________________________________________________________________
LogFile::~LogFile()
{

	fLog.close();

}


