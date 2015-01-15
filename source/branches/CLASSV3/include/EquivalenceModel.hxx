#ifndef _EQUIVALENCEMODEL_HXX
#define _EQUIVALENCEMODEL_HXX


/*!
 \file
 \brief Header file for EquivalenceModel class.
 
 
 @author BLG,BaM
 @version 2.0
 */

#include "IsotopicVector.hxx"
#include <math.h>
#include "CLASSObject.hxx"


using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a EquivalenceModel.
 The aim of these class is synthetyse all the commum properties to all 
 Equivalence Model.

 !!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!
Never instantiate EquivalenceModel in your CLASS input but it's derivated class
see @../Model/Equivalence/EQM_LIN_PWR_MOX.hxx  
see @../Model/Equivalence/EQM_QUAD_PWR_MOX.hxx
see @../Model/Equivalence/EQM_MLP_PWR_MOX.hxx

 @author BLG,BaM
 @version 3.0
 */
//________________________________________________________________________


class EquivalenceModel : public CLASSObject
{
	public :

	EquivalenceModel();
	EquivalenceModel(CLASSLogger* log);

	virtual ~EquivalenceModel();

	/// virtual method called to build a reprocessed fuel as a function of the burnup requierement the stock, mass....
	/*!
	 Build the fuel following the equivalance model with the proper requierment in term of mass, burnup....
	 \param double BurnUp desireted burnup reached by the fuel at the end of irradiation
	 \param double HMMass, needed Heavy metal mass needed
	 \param vector<double> &lambda, fraction of the stock to take (initialy should be 0)
	 \param vector<IsotopicVector> FissilArray, isotopicvectors to use to get the fissil part of the fuel
	 \param vector<IsotopicVector> FertilArray, isotopicvectors to use to get the fertil part of the fuel (if empty take it from the OutIncome)
	 */
	
	virtual	 vector<double> BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray );
	//}
	virtual void GuessLambda(vector<double>& lambda,int FirstStockID, int LastStockID, double DeltaM, vector<IsotopicVector> Stocks, double  HMMass);
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp) = 0; /*{return 0;}*/ //!< Return the molar fraction of fissile element in the fuel accodring to the Burnup, and a given fuel composition (this is the heart of the equivalence model) 

	

	IsotopicVector GetFertileList() {return fFertileList;}	//!<return the fertile list
	IsotopicVector GetFissileList() {return fFissileList;}	//!<return the fissile list

	void SetFertileList(IsotopicVector IV) {fFertileList = IV;}//!<set the fertile list
	void SetFissileList(IsotopicVector IV) {fFissileList = IV;}//!<set the fissile list
	/// Check either the IsotopicVector IV is in the validity domain of the models.
	/*!
	 return true if IV is in ValidityDomain
	 return false + a warning if IV is not in ValidityDomain
	 \param vector<IsotopicVector> IV, Fresh fuel composition
	 */
	//virtual  bool isIVInDomain(IsotopicVector IVFiss, double BU = 0 ) =0;
	
	protected :

	IsotopicVector fFertileList;	//!< contain the list of zai, needed as fertile, taken in a stock before fabrication
					//!< if no stock are provided, will take the isotopic vector in the Park income

	IsotopicVector fFissileList;	//!< contain the list of zai, needed as fissile, taken in a stock before fabrication
					//!< if no stock are provided the fuel will not be made


	private :

	void SetLambda(vector<double>& lambda ,int FirstStockID, int LastStockID, double LAMBDA_TOT);	//!< Set individual lambda according to the LAMBDA_TOT (lambda of all stocks)
	double FindLambdaMax( vector<IsotopicVector> Stocks, double  HMMass);							//!< Find the maximum LAMBDA_TOT of Stocks (ie lambda to reach HMass)
	double fOld_Lambda_Tot;//!<The old (old iteration) guessed lambda_tot (guessed from GuessLambda)
	double fLambda_max;//!<Value calculated by FindLambdaMax

};

#endif

