#ifndef _EQM_LIN_PWR_MOX_HXX
#define _EQM_LIN_PWR_MOX_HXX

#include "EquivalenceModel.hxx"

#include <string>

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
 based on a Linear Eq Model BU = SUM (alpha_i*n_i}

 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_LIN_PWR_MOX : public EquivalenceModel
{
	public :

	EQM_LIN_PWR_MOX(string WeightPath);
	~EQM_LIN_PWR_MOX();

	vector<double> BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray );

	private :

	string fWeightPath;
	vector<double> fFuelParameter;

};

#endif

