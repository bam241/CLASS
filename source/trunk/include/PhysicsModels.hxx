
#ifndef _PhysicsModels_HXX
#define _PhysicsModels_HXX


/*!
 \file
 \brief Header file for XS_INTERPOLATOR class.
 
 
 @authors BLG,BaM
 @version 1.0
 */
#include "EquivalenceModel.hxx"
#include "XSModel.hxx"
#include "IrradiationModel.hxx"
#include "EvolutionData.hxx"


using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
/*!
 Define all the Physics models used for a specific database. 
 
These class aim is basicaly to store 3 differents physics model : 
	
 The 2 following are data base related (for one Reactor and one fuel type ) :
User can either define his own (see manual) or uses the provided ones  :
\li XSModel : Mean cross section prediction (Closest, MLP )
\li EquivalenceModel : Fissile content prediction ( Linear,Quadratique, MLP , Baker & Ross, ...)

 This one is bateman solvers related :
\li IrradiationModel : can be Runge Kutta 4 or Matrix


 @authors BLG,BaM
 @version 1.0
 */
//________________________________________________________________________


class PhysicsModels : public CLASSObject
{

	public : 

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	PhysicsModels(); //!< Default Constructor
	
	//{
	/// XS, EM, IM Contructor
	/*!
	 Perform the Irradiation Calcultion using the XSSet data
	 \param XS : The XSModel (Mean cross section predictor)
	 \param EM : The EquivalenceModel (Fissile content predictor)
	 \param IM : The IrradiationModel (Bateman solver)
	 */
	PhysicsModels(XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM );
	//}
	
	//{
	/// CLASSLogger Contructor
	/*!
	 Perform the Irradiation Calcultion using the XSSet data
	 \param log : The CLASSLogger 
	 \param XS : The XSModel (Mean cross section predictor)
	 \param EM : The EquivalenceModel (Fissile content predictor)
	 \param IM : The IrradiationModel (Bateman solver)
	 */
	PhysicsModels(CLASSLogger* log, XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM );
	//}

	~PhysicsModels() {;}
	//{

	//{
	/// GenerateEvolutionData
	/*!
	 Call the 3 Physics models to compute the depletion calculation for the right fresh fuel with
	 the right mean cross sections
	 \param IV : The fresh fuel composition
	 \param cycletime : The irradiation time [s]
	 \param Power : The thermal (as always in CLASS) Power [W]
	 */
	EvolutionData GenerateEvolutionData(IsotopicVector IV, double cycletime, double Power);
	//}
	
	XSModel*		GetXSModel()   {return fXSModel;} //!< return the mean cross section predictor
	EquivalenceModel*	GetEquivalenceModel() {return fEquivalenceModel;} //!< return Fissile content predictor
	IrradiationModel*	GetIrradiationModel()  {return fIrradiationModel;} //!< return the Bateman solver
	
	PhysicsModels*		GetPhysicsModels()	{return this;}//!< return the PhysicsModels






 private :

 	XSModel* 		fXSModel;		//!< The XSModel (Mean cross sections prediction)
	EquivalenceModel*	fEquivalenceModel;	//!< The EquivalenceModel (Fresh fissile content prediction)
	IrradiationModel*	fIrradiationModel;	//!< The IrradiationModel (The Bateman's solver)



};

#endif

