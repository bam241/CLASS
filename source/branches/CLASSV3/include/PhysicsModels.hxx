
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
 Define all the physic models used for a specific database 
	
	The 2 following are data base related (for one Reactor and one fuel type and one ...) :
	User can either define his own (see manual) or uses the provided ones ) :

	XS_Interpolator = Closest , MLP  ....
	Equivalence_Model = Linear,Quadratique, MLP ... 

	this one is bateman solvers related : 
	(or it may be link to a evolution code (like MURE,DRAGON), not yet
	implemented but envisaged)
	IrradiationModel = RK4 or Matrix


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

	PhysicsModels();
	PhysicsModels(XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM );
	PhysicsModels(CLASSLogger* log, XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM );

	~PhysicsModels() {;}
	//{

	EvolutionData GenerateEvolutionData(IsotopicVector IV, double cycletime, double Power);

	XSModel*		GetXSModel()   {return fXSModel;}
	EquivalenceModel*	GetEquivalenceModel() {return fEquivalenceModel;}
	IrradiationModel*	GetIrradiationModel()  {return fIrradiationModel;}
	
	PhysicsModels*		GetPhysicsModels()	{return this;}






 private :

 	XSModel* 		fXSModel;              //!< The XSModel (Mean cross sections prediction)
	EquivalenceModel*	fEquivalenceModel; //!< The EquivalenceModel (Fresh fissile content prediction)
	IrradiationModel*	fIrradiationModel; //!< The IrradiationModel (The Bateman's solver)



};

#endif

