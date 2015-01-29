
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
#include "PhysicsModels.hxx"
#include "EvolutionData.hxx"

using namespace std;

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
	CLASSFuel(EvolutionData* evo);
	CLASSFuel(PhysicsModels* evo);


	virtual CLASSFuel* Clone()	{ return new CLASSFuel(*this); } //!< Correct way to copy a CLASSFuel in case of derivation


	EvolutionData* GetEvolutionData() {return fEvolutionData;}
	PhysicsModels* GetPhysicsModels() {return fPhysicsModels;}
	
	using CLASSObject::SetName;
	using CLASSObject::GetName;
	protected :


	private :

	EvolutionData* fEvolutionData;
	PhysicsModels* fPhysicsModels;
};

#endif

