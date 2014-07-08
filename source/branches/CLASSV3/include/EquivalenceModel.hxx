#ifndef _EQUIVALENCEMODEL_HXX
#define _EQUIVALENCEMODEL_HXX


/*!
 \file
 \brief Header file for EquivalenceModel class.
 
 
 @author BaM
 @version 2.0
 */

#include "IsotopicVector.hxx"
#include <math.h>
#include "CLASSObject.hxx"


using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a EquivalenceModel.
 The aim of these class is synthetyse all the commum properties to all Equivalence Model.
 
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________


class EquivalenceModel : public CLASSObject
{
	public :

	EquivalenceModel();
	EquivalenceModel(CLASSLogger* log);

	virtual ~EquivalenceModel();

	/// virtueal method called to build a reprocessed fuel as a function of the burnup requierement the stock, mass....
	/*!
	 Build the fuel following the equivalance model with the proper requierment in term of mass burnup....
	 \param double BurnUp desireted burnup reached by the fuel at the end of irradiation
	 \param double HMMass, needed Heavy metal mass needed
	 \param vector<double> &lambda, fraction of the stock to take (initialy should be 0)
	 \param vector<IsotopicVector> FissilArray, isotopicvectors to use to get the fissil part of the fuel
	 \param vector<IsotopicVector> FertilArray, isotopicvectors to use to get the fertil part of the fuel (if empty take it from the god)
	 */
	
	virtual	 vector<double> BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray );
	//}
	virtual void GuessLambda(vector<double>& lambda,int FirstStockID, int LastStockID, double DeltaM, vector<IsotopicVector> Stocks, double  HMMass);
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp) {return 0;}

	

	IsotopicVector GetFertileList() {return fFertileList;}
	IsotopicVector GetFissileList() {return fFissileList;}

	void SetFertileList(IsotopicVector IV) {fFertileList = IV;}
	void SetFissileList(IsotopicVector IV) {fFissileList = IV;}

	
	protected :

	IsotopicVector fFertileList;	//!< contain the list of zai, needed as fertile, taken in a stock before fabrication
					//!< if no stock are provided, will take the isotopic vector in the Park income

	IsotopicVector fFissileList;	//!< contain the list of zai, needed as fissile, taken in a stock before fabrication
					//!< if no stock are provided, will take the isotopic vector in the Park income


	private :

	void SetLambda(vector<double>& lambda ,int FirstStockID, int LastStockID, double LAMBDA_TOT);
	void FindLambdaMax(int FirstStockID, int LastStockID, vector<IsotopicVector> Stocks, double  HMMass);
	double fOld_Lambda_Tot[2];
	double fLambda_max[2];

};

#endif

