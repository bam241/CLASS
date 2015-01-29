
#ifndef _CLASSOBJECT_HXX
#define _CLASSOBJECT_HXX


/*!
 \file
 \brief Header file for CLASSObject class.
 
 
 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include "CLASSLogger.hxx"

#include "TNamed.h"

using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a CLASS Object.
 The aim of these class is synthetyse all the commum properties to all CLASS Element.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class CLASSObject : public TNamed
{
public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	//{
	/// Normal Constructor.
	/*!
	 Make a new CLASSObject
	 */
	CLASSObject();
	//}

	//{
	/// Log Constructor.
	/*!
	 Make a new CLASSObject
	 \param log : used for the log.
	 */
	
	CLASSObject(CLASSLogger* log);
	//}
	//@}

	/*!
	 \name Clone
	 */
	//@{

	virtual CLASSObject* Clone()	{ return new CLASSObject(*this); } //!< Correct way to copy a CLASSObject in case of derivation
	//}
	//@}

	
	/*!
	 \name Set/Get
	 */
	//@{


#ifndef __CINT__
	void		SetLog(CLASSLogger* log)	{ fLog = log;}		//!< Set the CLASSLogger
	CLASSLogger*	GetLog()		{ return fLog; }		//!< Return the Pointer to the Log
#endif
	//@}

	using TNamed::SetName;
	using TNamed::GetName;
protected :
#ifndef __CINT__
	CLASSLogger*	fLog;			//!< Pointer to the Log
#endif


private :

	ClassDef(CLASSObject,0);
};

#endif
