
#ifndef _LOGFILE_CLASS
#define _LOGFILE_CLASS


/*!
 \file 
 \brief Header file for LogFile class. 

 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>
#include "stdlib.h"
using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a LogFile.
 The aim of this class is to centralize the all CLASS software message inside a file.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class LogFile
{
public:

	//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{


	LogFile(string LogFileName );	//!< Normal Constructor

	~LogFile();	//!< Normal Destructor

	//@}

	//********* In/Out Method *********//

	/*!
	 \name In/Out
	 */
	//@{
	string GetLogFileName() const { return fLogFileName; }	//!w return the logfile name

	std::ofstream fLog;		//!< Log Stream

	//@}

	private :

	string fLogFileName;		//!< Log File name
};

#endif



