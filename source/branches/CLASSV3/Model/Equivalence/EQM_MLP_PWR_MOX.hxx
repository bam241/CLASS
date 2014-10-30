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

	EQM_MLP_MOX(string TMVAWeightPath);	//!< Constructor  string TMVAWeightPath => PATH/TMVAWeight.xml (path to tmva weight)
	EQM_MLP_MOX(CLASSLogger* log, string TMVAWeightPath);	//!< Constructor CLASSLogger* log ,string TMVAWeightPath => PATH/TMVAWeight.xml

	double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp); //!<Return the molar fraction of fissile element thanks to a Multi Layer Perceptron

	private :
	TTree* CreateTMVAInputTree(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);//!<Create input tmva tree to be read by ExecuteTMVA
	double ExecuteTMVA(TTree* theTree);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree


	string fTMVAWeightPath;;//!<The weight needed by TMVA to construct and execute the multilayer perceptron

};

#endif

