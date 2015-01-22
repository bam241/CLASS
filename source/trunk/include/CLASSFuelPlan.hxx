
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
 The aim of these class is synthetyse all the commum properties to all CLASS Fuel Element.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class CLASSFuelPlan : public CLASSObject
{
	public :
	///< Normal Constructor.
	CLASSFuelPlan();
	CLASSFuelPlan(CLASSLogger* log);

	void AddFuel(cSecond time,  CLASSFuel fuel, double BurnUp);
	
	void AddFuel(cSecond time,  EvolutionData* fuel, double BurnUp) {AddFuel( time, CLASSFuel(fuel), BurnUp);}
	void AddFuel(cSecond time,  PhysicsModels* fuel, double BurnUp) {AddFuel( time, CLASSFuel(fuel), BurnUp);}

	pair< CLASSFuel, double> GetFuelAt(cSecond t);

	using CLASSObject::SetName;
	using CLASSObject::GetName;
	protected :


	private :

	map< cSecond, pair< CLASSFuel, double > >	fLoadingPlan;	///< Loading PLan to change the EvolutionData (and the associetedBurnup) according to the Plan

};

#endif

