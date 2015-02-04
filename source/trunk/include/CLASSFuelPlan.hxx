
#ifndef _CLASSFUELPLAN_HXX
#define _CLASSFUELPLAN_HXX


/*!
 \file
 \brief Header file for CLASSFuelPlan class.


 @author BaM
 @version 2.0
 */

#include <string>
#include <fstream>
#include <map>

#include "CLASSObject.hxx"
#include "CLASSFuel.hxx"

using namespace std;
typedef long long int cSecond;


//-----------------------------------------------------------------------------//
/*!
 Define a CLASS Object.
 The aim of these class is to allow a Reactor to change its CLASSFuel and/or burnup
 during its lifetime

 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class CLASSFuelPlan : public CLASSObject
{
	public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	//{
	/// CLASSFuelPlan Constructor.
	/*!
	 Make a new CLASSFuelPlan
	 */
	CLASSFuelPlan();
	//}
	
	//{
	/// CLASSFuelPlan Constructor.
	/*!
	 Make a new CLASSFuelPlan
	 \param log : used for the log.
	 */	CLASSFuelPlan(CLASSLogger* log);
	//}
	//@}

	/*!
	 \name Adding Method
	 */
	//@{
	
	void AddFuel(cSecond time,  CLASSFuel fuel, double BurnUp);	//!< Add A new CLASSFuel at the corresponding time and Burnup
	void AddFuel(cSecond time,  EvolutionData* fuel, double BurnUp)
			{ AddFuel( time, CLASSFuel(fuel), BurnUp); }	//!< Add A new EvolutionData at the corresponding time and Burnup
	void AddFuel(cSecond time,  PhysicsModels* fuel, double BurnUp)
			{ AddFuel( time, CLASSFuel(fuel), BurnUp); }	//!< Add A new Physicis Model at the corresponding time and Burnup

	//a}
	
	/*!
	 \name Get Method
	 */
	//@{

	pair< CLASSFuel, double> GetFuelAt(cSecond t);			//!< Get fuel and associated Burnup loaded at the time t

	
	//@}
	using CLASSObject::SetName;
	using CLASSObject::GetName;
	protected :


	private :

	map< cSecond, pair< CLASSFuel, double > >	fLoadingPlan;	///< Loading plan to change the EvolutionData (and the associeted Burnup) according to the plan

};

#endif

