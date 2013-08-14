
#ifndef _LOGFILE_CLASS
#define _LOGFILE_CLASS


/*!
 \file 
 \brief Header file for LogFile class. 
  The aim of this Class is to centralize the all CLASS software message inside a file.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>
#include "stdlib.h"
using namespace std;


class LogFile
{
public:
	//!< Normal Constructor
	LogFile(string LogFileName );

	//!< Normal Destructor
	~LogFile();
	
	string GetLogFileName() const { return fLogFileName; }

	std::ofstream fLog;		//!< Log Stream
	
	private :
	string fLogFileName;
};

#endif

