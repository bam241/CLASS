#ifndef _ZAIMass_
#define _ZAIMass_

/*!
 \file
 \brief Header file for ZAIMass classes.
 
 
 @author BaM
 @version 2.0
 */

#include <map>
#include "ZAI.hxx"
#include "TObject.h"
#include <iostream>

using namespace std;




///< A ZAIMass .

class ZAIMass
{
	
	
public:
		///< Default constructor
	ZAIMass();
		///< Normal Constructor.

		///< Normal Destructor.
	~ZAIMass();
	
	map<ZAI, double> fZAIMass; //! ZAI mass list


	double GetMass(const int Z, const int A ) const;

};


#endif
