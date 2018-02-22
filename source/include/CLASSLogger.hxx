#ifndef _LOG_
#define _LOG_


/*!
 \file 
 \brief Header file for CLASSLogger class. 

 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include <iostream>
#include <cstring>
#include <sstream>
#include "stdlib.h"
using namespace std;


#ifndef __ROOTCLING__

#define ERROR		if(fLog)if(fLog->GetMaxOutPutLVL() >= 0) fLog->E() << "!!!ERROR!!! " << "[" << __FILE__ << ":" << __FUNCTION__ << "] "
#define WARNING		if(fLog)if(fLog->GetMaxOutPutLVL() >= 1) fLog->W() << "!!WARNING!! " << "[" << __FILE__ << ":" << __FUNCTION__ << "] "
#define INFO		if(fLog)if(fLog->GetMaxOutPutLVL() >= 2) fLog->I() << "!!!!INFO!!! " << "[" << __FILE__ << "] "

#define DBGL		if(fLog)if(fLog->GetMaxOutPutLVL() >= 3) fLog->D() << __FILE__ << " : " << __LINE__ << " [" << __FUNCTION__ << "]" << endl;
#define DBGV(x)		{if(fLog)if(fLog->GetMaxOutPutLVL() >= 3) fLog->D() << __FILE__ << " : " << __LINE__ << " [" << __FUNCTION__ << "]" << x << endl;}

#else

#define ERROR		cout
#define INFO		cout
#define WARNING		cout
#define DBGL
#define DBGV(x)


#endif



#ifndef _LOGTYPE_
#define _LOGTYPE_

//-----------------------------------------------------------------------------//
//!handles output stream in CLASS

/*!
 Define a LogType.
 The aim of this class is to handle output stream in CLASS.
 
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________


class LogType
{
public:
	//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	//{
	/// Normal Constructor.
	/*!
	 Make a new LogType
	 \param Log : ostream output
	 */
	LogType(ostream &Log) { fLog = &Log; fLog2 = 0; }	//!< Normal Constructor
	//}
	
	~LogType()  {}	//!< Normal Destructor

	//@}

	//********* In/Out Method *********//

	/*!
	 \name In/Out
	 */
	//@{
	string GetCLASSLoggerName() const { return fCLASSLoggerName; }	//!< return the CLASSLogger name

	LogType &operator << (std::ostream& (*manip)(std::ostream &))
	{
		manip( *(this->fLog) );
		if(fLog2)
			manip( *(this->fLog2) );
		return *this;
	}


	template<typename T>
	inline LogType& operator << (T something)
	{
                *(this->fLog) << something;
		if(fLog2)
			*(this->fLog2) << something;
		 return *this;
	}
	//}

	void SetSecondOutput(ostream &log) {fLog2 = &log;} // used to direct the stream into two output ostream

	private :

	ostream *fLog;				//!< Main ostream output
	ostream *fLog2;				//!< secondary ostream output

	string fCLASSLoggerName;		//!< Log File name
};


#endif


#ifndef _CLASSLogger_
#define _CLASSLogger_

//-----------------------------------------------------------------------------//
//! Object to handle output messages

/*!
 Define a CLASSLogger.
 The aim of this class is to centralize the all CLASS software message inside a file.
 
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________


class CLASSLogger
{
public:

	//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	//{
	/// Normal Constructor.
	/*!
	 Make a new LogType
	 \param CLASSLoggerName : name of the CLASSLogFile wnated for the log
	 \param VerboseLvl : verbose level in terminal
	 \param OutputLvl : verbose level in the CLASSLogFile
	 
	 The different available levels are :
	 \li 0 : ERROR only
	 \li 1 : WARNING + lvl 0
	 \li 2 : INFO + lvl 1
	 \li 3 : DEBUG + lvl 2
	 
	 */
	CLASSLogger(string CLASSLoggerName = "CLASS_OUTPUT.log", int VerboseLvl = 0, int OutputLvl = 1 );
	//}
	
	~CLASSLogger();	//!< Normal Destructor

	//@}

	//********* In/Out Method *********//

	/*!
	 \name In/Out
	 */
	//@{
	string GetCLASSLoggerName() const { return fCLASSLoggerName; }	//!< return the CLASSLogger name
	int GetMaxOutPutLVL()	const { return fMaxOutPutLVL; }		//!< Return File Output lvl
	int GetVerboseLVL()	const { return fVerboseLVL; }

	LogType E() {return *fError;}		//!< Return the ERROR Streamer
	LogType W() {return *fWarning;}		//!< Return the WARNING Streamer
	LogType I() {return *fInfo;}		//!< Return the INFO Streamer
	LogType D() {return *fDebug;}		//!< Return the DEBUG Streamer

//@}



	
	private :
	int fMaxOutPutLVL;			//!< Maximal output/verbose lvl
	int fVerboseLVL;

	LogType* fError;			//!< ERROR streamer
	LogType* fInfo;				//!< INFO streamer
	LogType* fWarning;			//!< WARNING streamer
	LogType* fDebug;			//!< DEBUG streamer

	ofstream fOutPutFile;			//!< Log Output File name
	string fCLASSLoggerName;		//!< Log File name
};

#endif



#endif



