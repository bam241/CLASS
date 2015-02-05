
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
//!  Allows to define PhysicsModels & EvolutionData as a CLASSFuel

/*!
 Define a CLASS Object.
 The aim of these class is to handle PhysicsModels &
 EvolutionData as a CLASSFuel .
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________



class CLASSFuel : public CLASSObject
{
	public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	//{
	/// EvolutionData Constructor.
	/*!
	 Make a new CLASSObject
	 /param evo : EvolutionData stored
	 */
	CLASSFuel(EvolutionData* evo);
	//}
	
	//{
	/// PhysicsModels Constructor.
	/*!
	 Make a new CLASSObject
	 /param evo : PhysicsModels stored
	 */
	CLASSFuel(PhysicsModels* evo);
	//}
	//@}
	
	
	virtual CLASSFuel* Clone()	{ return new CLASSFuel(*this); } //!< Correct way to copy a CLASSFuel in case of derivation
	
	
	EvolutionData* GetEvolutionData() {return fEvolutionData;}	//!< Return the EvolutionData (NULL if PhysicsModels)
	PhysicsModels* GetPhysicsModels() {return fPhysicsModels;}	//!< Return the PhysicsModels (NULL if EvolutionData)
	
	using CLASSObject::SetName;
	using CLASSObject::GetName;
	protected :
	
	
	private :
	
	EvolutionData* fEvolutionData;		//!< the EvolutionData (NULL if PhysicsModels)
	PhysicsModels* fPhysicsModels;		//!< the PhysicsModels (NULL if EvolutionData)
};

#endif

