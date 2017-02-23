
#ifndef _REACTORSCHEDULER_HXX
#define _REACTORSCHEDULER_HXX


/*!
 \file
 \brief Header file for ReactorScheduler class.


 @author BLG, BaM
 @version 2.0
 */

#include <string>
#include <fstream>
#include <map>
#include "EvolutionData.hxx"
#include "PhysicsModels.hxx"

#include "CLASSObject.hxx"

using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
//!  Allows to define PhysicsModels & EvolutionData as a ReactorModel

/*!
 Define a CLASS Object.
 The aim of these class is to handle PhysicsModels &
 EvolutionData as a ReactorModel .
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________
class ReactorModel : public CLASSObject
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
	ReactorModel(EvolutionData* evo);
	//}
	
	//{
	/// PhysicsModels Constructor.
	/*!
	 Make a new CLASSObject
	 /param evo : PhysicsModels stored
	 */
	ReactorModel(PhysicsModels* evo);
	//}
	//@}
	
	
	virtual ReactorModel* Clone()	{ return new ReactorModel(*this); } //!< Correct way to copy a ReactorModel in case of derivation
	
	
	EvolutionData* GetEvolutionData() {return fEvolutionData;}	///< Return the EvolutionData (NULL if PhysicsModels)
	PhysicsModels* GetPhysicsModels() {return fPhysicsModels;}	///< Return the PhysicsModels (NULL if EvolutionData)
	
	using CLASSObject::SetName;
	using CLASSObject::GetName;

	protected :
	
	private :
	
	EvolutionData* fEvolutionData;		///< the EvolutionData (NULL if PhysicsModels)
	PhysicsModels* fPhysicsModels;		///< the PhysicsModels (NULL if EvolutionData)

};
//----------------------------------------------------------------------------//
/*!
 Define a CLASS Object.
 The aim of these class is to define a SchedulerEntry. This object contains 
 informations (ReactorModel, BurnUpPower, HMMass) of one entry of the 
 ReactorScheduler (@see ReactorScheduler)

 @author  BLG
 @version 1.0
 */
//________________________________________________________________________
class ScheduleEntry : public CLASSObject
{
  public:
  	/*!
	 \name Constructor/Desctructor
	 */
	//@{
  	ScheduleEntry(ReactorModel* ED_Or_PhysMod, double BurnUp, double Power, double HMMass );
  	//@}

  	/*!
	 \name Get methods
	 */
	//@{
    ReactorModel* GetReactorModel(){return fReactorModel;}	///< Return the Reactor Model
    double GetBurnUp(){return fBurnUp;}	///< Return the Reactor BurnUp
    double GetPower(){return fPower;}	///< Return the Reactor Power
    double GetHeavyMetalMass(){return fHeavyMetalMass;} ///< Return the Reactor Heavy Metal Mass
	//@}

  private:

  	ReactorModel* fReactorModel ; 	///< the Reactor Model
    double fBurnUp; 	///< the Reactor BurnUp
    double fPower;       	///< the Reactor Power  
    double fHeavyMetalMass;	///< the Reactor Heavy Metal Mass

};

//----------------------------------------------------------------------------//
//! Allows a Reactor to change its ED_Or_PhysMod,  BurnUp,  Power,  HMMass

/*!
 Define a CLASS Object.
 The aim of these class is to allow a Reactor to change its  ED_Or_PhysMod,  
 BurnUp,  Power,  HMMass during its lifetime

 @author BaM, BLG
 @version 2.0
 */
//________________________________________________________________________

class ReactorScheduler : public CLASSObject
{
	public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	//{
	/// ReactorScheduler Constructor.
	/*!
	 Make a new ReactorScheduler
	 */
	ReactorScheduler();
	//}
	
	//{
	/// ReactorScheduler Constructor.
	/*!
	 Make a new ReactorScheduler
	 \param log : used for the log.
	 */	
	 ReactorScheduler(CLASSLogger* log);
	//}
	//@}

	/*!
	 \name Adding Method
	 */
	//@{
	
	void AddEntry(cSecond time,  ReactorModel*   Model, double BurnUp, double Power, double HMMass);	//!< Add A new Entry to the scheduler 

	//@}
	
	/*!
	 \name Get Method
	 */
	//@{
	ScheduleEntry* GetEntryAt(cSecond t); ///< Get Schedule entry at time t
	//@}

	using CLASSObject::SetName;
	using CLASSObject::GetName;

	protected :


	private :

	map< cSecond, ScheduleEntry* >	fReactorSchedulerMap;	///< Get the reactor scheduler map

};

#endif

