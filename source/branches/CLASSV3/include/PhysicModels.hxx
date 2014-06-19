
#ifndef _PHYSICMODELS_HXX
#define _PHYSICMODELS_HXX


/*!
 \file
 \brief Header file for XS_INTERPOLATOR class.
 
 
 @authors BLG,BaM
 @version 1.0
 */
 #include "CLASSObject.hxx"
 #include "EquivalenceModel.hxx"
 #include "XSModel.hxx"
 #include "IrradiationModel.hxx"


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
	(or it may be link to a evolution code (like MURE), not yet implemented but envisaged)
	IrradiationModel = RK4 or Matrix


 @authors BLG,BaM
 @version 1.0
 */
//________________________________________________________________________


class PhysicModels : public CLASSObject
{

	public : 

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	PhysicModels(XSModel XS,EquivalenceModel EM,IrradiationModel IM );	

	~PhysicModels();
	//{

	EvolutionData GenerateEvolutionData(IsotopicVector IV, double cycletime, double Power);

	XSModel GetXSModel()   {return fXSModel;}
	EquivalenceModel GetEquivalenceModel() {return fEquivalenceModel;}
	IrradiationModel GetIrradiationModel()  {return fIrradiationModel;}






 private :

 	XSModel 			fXSModel;
	EquivalenceModel	fEquivalenceModel;
	IrradiationModel	fIrradiationModel;



};

#endif

