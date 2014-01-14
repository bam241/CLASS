
#ifndef _CLASSOBJECT_HXX
#define _CLASSOBJECT_HXX


/*!
 \file
 \brief Header file for CLSSObject class.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include "LogFile.hxx"

#include "TNamed.h"

using namespace std;


class CLSSObject : public TNamed
{
public :
	///< Normal Constructor.
	CLSSObject();
	virtual CLSSObject* Clone()	{ return new CLSSObject(*this); } //!< Correct way to copy a CLSSObject in case of derivation


	void		SetLog(LogFile* log)	{ fLog = log; }
	void		IsLog(bool islog)	{ fNoLog = islog; }
	
	LogFile*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log
	bool		PrintLog()			{ return fNoLog; }
	
private :
 	LogFile*	fLog;			//!< Pointer to the Log
	bool		fNoLog;
	ClassDef(CLSSObject,0);
};

#endif

