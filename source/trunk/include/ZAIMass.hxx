#ifndef _ZAIMass_
#define _ZAIMass_

/*!
 \file
 \brief Header file for ZAIMass classes.
 
 
 @author BaM & BaL
 @version 2.0
 */

#include <map>

#include "ZAI.hxx"
#include "TObject.h"
#include <iostream>

using namespace std;


class IsotopicVector;

//-----------------------------------------------------------------------------//
//! Defines the molar mass of a ZAI

/*!
 The aims of this class is to handle the molar mass of each ZAI
 
 @author BaM, BaL
 @version 1.0
 */
//________________________________________________________________________


class ZAIMass
{
	
	
public:
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	ZAIMass();//!< Normal Constructor.
	
	~ZAIMass();//!< Normal Destructor.
	//@}

	/*!
	 \name Fucntions returning molar mass [g/mol]
	 */
	//@{
	double GetMass(ZAI zai ) const; //!< get with ZAI
	double GetMass(const int Z, const int A )    const { return GetMass( ZAI(Z, A, 0) ); } //!< Get with Z, A
	//@}
	
	double GetMass(const IsotopicVector IV)    const; //return Mass of IV [t]

private:
	map<ZAI, double> fZAIMass; //! ZAI mass list


};


#endif
