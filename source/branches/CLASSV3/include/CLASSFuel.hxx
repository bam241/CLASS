
#ifndef _CLASSFUEL_HXX
#define _CLASSFUEL_HXX


/*!
 \file
 \brief Header file for CLASSFuel class.


 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>

#include "CLASSObject.hxx"

using namespace std;

class EvolutionData;
class PhysicsModels;
//-----------------------------------------------------------------------------//
/*!
 Define a CLASS Object.
 The aim of these class is synthetyse all the commum properties to all CLASS Fuel Element.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class CLASSFuel : public CLASSObject
{
	public :
	///< Normal Constructor.
	CLASSFuel();
	CLASSFuel(CLASSLogger* log);


	virtual CLASSFuel* Clone()	{ return new CLASSFuel(*this); } //!< Correct way to copy a CLASSFuel in case of derivation


	virtual EvolutionData* GetEvolutionData() {return 0;}
	virtual PhysicsModels* GetPhysicsModels() {return 0;}
	using CLASSObject::SetName;
	using CLASSObject::GetName;
	protected :


	private :

};

#endif

