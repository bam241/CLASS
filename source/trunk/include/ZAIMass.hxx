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

///< A ZAIMass .

class ZAIMass
{
	
	
public:
		///< Default constructor
	ZAIMass();
		///< Normal Constructor.

		///< Normal Destructor.
	~ZAIMass();
	


	double GetMass(ZAI zai ) const;
	double GetMass(const int Z, const int A )    const { return GetMass( ZAI(Z, A, 0) ); }

	double GetMass(const IsotopicVector IV)    const; //return Mass of IV in tons

private:
	map<ZAI, double> fZAIMass; //! ZAI mass list


};


#endif
