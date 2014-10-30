#ifndef _LOG_CLASS_
#define _LOG_CLASS_


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


#ifndef __CINT__

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



//-----------------------------------------------------------------------------//
/*!
 Define a CLASSLogger.
 The aim of this class is to centralize the all CLASS software message inside a file.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________
#ifndef _LOGTYPE_CLASS
#define _LOGTYPE_CLASS



class LogType
{
public:
	//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{


	LogType(ostream &Log) { fLog = &Log; fLog2 = 0; }	//!< Normal Constructor

	~LogType()  {}	//!< Normal Destructor

	//@}

	//********* In/Out Method *********//

	/*!
	 \name In/Out
	 */
	//@{
	string GetCLASSLoggerName() const { return fCLASSLoggerName; }	//!w return the CLASSLogger name

	LogType &operator<<(std::ostream& (*manip)(std::ostream &))
	{
		manip( *(this->fLog) );
		if(fLog2)
			manip( *(this->fLog2) );
		return *this;
	}


	template<typename T>
	inline LogType& operator<<(T something)
	{
                *(this->fLog) << something;
		if(fLog2)
			*(this->fLog2) << something;
		 return *this;
	}


	void SetSecondOutput(ostream &log) {fLog2 = &log;}

	private :

	ostream *fLog;
	ostream *fLog2;

	string fCLASSLoggerName;		//!< Log File name
};


#endif


#ifndef _CLASSLogger_CLASS_
#define _CLASSLogger_CLASS_



class CLASSLogger
{
public:

	//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	CLASSLogger();

	CLASSLogger(string CLASSLoggerName, int VerboseLvl = 0, int OutputLvl = 1 );	//!< Normal Constructor

	~CLASSLogger();	//!< Normal Destructor

	//@}

	//********* In/Out Method *********//

	/*!
	 \name In/Out
	 */
	//@{
	string GetCLASSLoggerName() const { return fCLASSLoggerName; }	//!w return the CLASSLogger name
	int GetMaxOutPutLVL()	const { return fMaxOutPutLVL; }
	int GetVerboseLVL()	const { return fVerboseLVL; }

	LogType E() {return *fError;}
	LogType W() {return *fWarning;}
	LogType D() {return *fDebug;}
	LogType I() {return *fInfo;}

//@}



	
	private :
	int fMaxOutPutLVL;
	int fVerboseLVL;

	LogType* fError;
	LogType* fInfo;
	LogType* fWarning;
	LogType* fDebug;

	ofstream fOutPutFile;
	string fCLASSLoggerName;		//!< Log File name
};

#endif



#endif



