
#ifndef _LOGFILE_CLASS
#define _LOGFILE_CLASS


/*!
 \file 
 \brief Header file for LogFile class. 
*/
#include <string>
#include <fstream>
#include "stdlib.h"
using namespace std;


class LogFile
{
	public:
		LogFile(string LogFileName );	//!< Normal Constructor
		~LogFile();			//!< Normal Destructor

		std::ofstream fLog;		//!< Log Stream

};

#endif

