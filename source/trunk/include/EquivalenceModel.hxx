#ifndef _EQUIVALENCEMODEL_
#define _EQUIVALENCEMODEL_


/*!
 \file
 \brief Header file for EquivalenceModel class.
 
 
 @author BaM
 @author BLG
 @version 3.0
 */

#include "IsotopicVector.hxx"
#include <math.h>
#include "CLASSObject.hxx"


using namespace std;

//-----------------------------------------------------------------------------//

//! Determines how to build a fresh fuel 
/*!
 Define an EquivalenceModel.
 The aim of these class is to gather all the commum properties of all
 Equivalence Model.
 
 \warning
 Never instantiate EquivalenceModel in your CLASS input but it's derivated class
 @see EQM_BakerRoss_FBR_MOX
 @see EQM_LIN_PWR_MOX
 @see EQM_MLP_PWR_MOX
 @see EQM_POL_PWR_UO2
 @see EQM_QUAD_PWR_MOX
 
 @author BaM
 @author BLG
 @version 3.0
 */
//________________________________________________________________________


class EquivalenceModel : public CLASSObject
{
	public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	EquivalenceModel();			//!< Default constructor
	EquivalenceModel(CLASSLogger* log);	//!< Logger constructor
	
	virtual ~EquivalenceModel();		//!< Destructor
	//@}
	
	/*!
	 \name Fuel Construction Method
	 */
	//@{
	
	//{
	/// BuildFuel function.
	/*!
	 Build the fuel following the equivalance model with the proper requierment in term of mass, burnup....
	 \param double burnup reached by the fuel at the end of irradiation
	 \param double HMMass, Heavy metal mass needed
	 \param vector<double> &lambda, fraction of the stock to take (initialy should be 0)
	 \param vector<IsotopicVector> FissilArray, isotopicvectors to use to get the fissile part of the fuel
	 \param vector<IsotopicVector> FertilArray, isotopicvectors to use to get the fertile part of the fuel (if empty take it from the OutIncome)
	 */
	virtual	 vector<double> BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray );
	//}
	
	virtual void GuessLambda(vector<double>& lambda,
				 int FirstStockID, int LastStockID, double DeltaM,
				 vector<IsotopicVector> Stocks, double  HMMass); //!< Guess a portion of IsotopicVectors to take (dichotomy)
	
	//@}
	
	/*!
	 \name Get/Set Method
	 */
	//@{
	
	IsotopicVector GetFertileList() {return fFertileList;}	//!<return the fertile list
	IsotopicVector GetFissileList() {return fFissileList;}	//!<return the fissile list
	
	double GetBuildFuelFirstGuess(){return fFirstGuessFissilContent;} //!<Get the initialization value for BuildFuel algorithm
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp) {return 0;} //!< Return the molar fraction of fissile element in the fuel according to the burnup, and a given fuel composition (this is the heart of the equivalence model)

	double GetRelativMassPrecision() const	{ return fRelativMassPrecision; }	//!< Mass precision
	int GetMaxInterration()		 const	{ return fMaxInterration; }		//!< Max iterration in build fueld algorythm

	double GetActualFissileContent() const { return fActualFissileContent; }	//!< Get the fissile content at the actual dichotomy step (usefull for EQM_MLP_Kinf)
	
	void SetFertileList(IsotopicVector IV) {fFertileList = IV;}//!<set the fertile list
	void SetFissileList(IsotopicVector IV) {fFissileList = IV;}//!<set the fissile list
	void SetBuildFuelFirstGuess(double FirstGuess){fFirstGuessFissilContent = FirstGuess;} //!<set the initialization value for BuildFuel algorithm
	void SetRelativMassPrecision( double val)	{ fRelativMassPrecision = val; }	//!< Mass precision
	void SetMaxInterration(int val)			{ fMaxInterration = val; }		//!< Max iterration in build fueld algorythm

	
	//@}
	
	
	
	
	
	protected :
	
	IsotopicVector fFertileList;	//!< contain the list of zai, needed as fertile, taken in a stock before fabrication
	//!< if no stock are provided, will take the isotopic vector in the Park income
	
	IsotopicVector fFissileList;	//!< contain the list of zai, needed as fissile, taken in a stock before fabrication
	//!< if no stock are provided the fuel will not be made
	
	double fFirstGuessFissilContent;//!< fissile content for BuildFuel initialization (in weight proportion)
	double fActualFissileContent;	//!< fissile content at the actual dichotomy step (usefull for EQM_MLP_Kinf)
	
	
	
	private :
	
	double fOld_Lambda_Tot_Minus;//!< The old (old iteration) guessed lambda_tot (guessed from GuessLambda)
	double fOld_Lambda_Tot_Plus ;//!< The old (old iteration) guessed lambda_tot (guessed from GuessLambda)
	double fLambda_max;		//!< Value calculated by FindLambdaMax
	
	double fRelativMassPrecision;	//!< Mass precision
	int fMaxInterration;		//!< Max iterration in build fueld algorythm

	

	double FindLambdaMax( vector<IsotopicVector> Stocks, double  HMMass);	//!< Find the maximum LAMBDA_TOT of Stocks (ie lambda to reach HMass)
	void SetLambda(vector<double>& lambda ,int FirstStockID, int LastStockID, double LAMBDA_TOT);	//!< Set individual lambda according to the LAMBDA_TOT (lambda of all stocks)


};

#endif









