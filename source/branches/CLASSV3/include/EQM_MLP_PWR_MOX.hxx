#ifndef _EQM_MLP_MOX_HXX
#define _EQM_MLP_MOX_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"

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


class EQM_MLP_MOX : public EquivalenceModel
{
	public :
	
	EQM_MLP_MOX(string TMVAWeightPath);
	EQM_MLP_MOX(CLASSLogger* log, string TMVAWeightPath);

	double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);

	private :
	TTree* CreateTMVAInputTree(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);
	double ExecuteTMVA(TTree* theTree);


	string fTMVAWeightPath;

};

#endif

