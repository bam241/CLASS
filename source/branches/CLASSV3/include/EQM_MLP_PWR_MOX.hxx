#ifndef _EQM_MLP_MOX_HXX
#define _EQM_MLP_MOX_HXX

#include "CLASSHeaders.hxx"

/*!
 \file
 \brief Header file for EQM_MLP_MOX class.
 
 
 @author BaM
 @version 1.0
 */


using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a EQM_MLP_MOX.
 The aim of these class is to constuct a fuel from an equivalence model
 based on a  Multi layer perceptron

 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class IsotopicVector;


class EQM_MLP_MOX  : public EquivalenceModel
{
	public :
	
	EQM_MLP_MOX(string TMVAWeightPath);

	///  method called to build a reprocessed fuel as a function of the burnup requierement the stock, mass....
	/*!
	 Build the fuel following the equivalance model with the proper requierment in term of mass burnup....
	 \param double BurnUp desireted burnup reached by the fuel at the end of irradiation
	 \param double HMMass, needed Heavy metal mass needed
	 \param vector<double> &lambda, fraction of the stock to take (initialy should be 0)
	 \param vector<IsotopicVector> FissilArray, isotopicvectors to use to get the fissil part of the fuel
	 \param vector<IsotopicVector> FertilArray, isotopicvectors to use to get the fertil part of the fuel (if empty take it from the god)
	 */
	
	vector<double>  BuildFuel(double BurnUp, double HMMass,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray = vector<IsotopicVector>());
	//}
	
	private :
	void CreateTMVAInputTree(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);
	double ExecuteTMVA();
	void GuessLambda(double& lambda, int& StockID,int FirstStockID, int LastStockID, double DeltaM,double StockHM);


	string fTMVAWeightPath;

};

#endif

