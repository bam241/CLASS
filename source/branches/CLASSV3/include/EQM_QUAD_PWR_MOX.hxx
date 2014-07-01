#ifndef _EQM_QUAD_PWR_MOX_HXX
#define _EQM_QUAD_PWR_MOX_HXX

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
 based on a Quadratic Pu equivalent Model

 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_QUAD_PWR_MOX : public EquivalenceModel
{
	public :
	
	EQM_QUAD_PWR_MOX(string WeightPath);
	~EQM_QUAD_PWR_MOX();

	double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);

	private :

	string fWeightPath;
	vector<double> fFuelParameter;

};

#endif

