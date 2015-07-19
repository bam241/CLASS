#ifndef _EQM_PWR_MLP_MOX_AM_HXX
#define _EQM_PWR_MLP_MOX_AM_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"

using namespace std;

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on neural network

/*!
 The aim of these class is to constuct a fuel from an equivalence model
 based on a  Multi layer perceptron

 @author BLG
 @version 3.0
 */
//________________________________________________________________________

class EQM_PWR_MLP_MOX_AM : public EquivalenceModel
{
	public :
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// normal constructor
	/*!
	 Create a EQM_PWR_MLP_MOX_AM 
	 \param  TMVAWeightPath :  PAth to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 */
	EQM_PWR_MLP_MOX_AM(string TMVAWeightPath);
	//}
	
	//{
	/// Logger constructor
	/*!
	 Create a EQM_PWR_MLP_MOX_AM
	 \param log : use for log
	 \param  TMVAWeightPath :  PAth to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 */
	EQM_PWR_MLP_MOX_AM(CLASSLogger* log, string TMVAWeightPath);
	//}
	//@}
	
	//{
	/// Return the molar fissile fraction according fissile & ferile content using a Multi Layer Peceptron (MLP)
	/*!
	 \param Fissil : The composition of the fissile matter
	 \param Fertil : The composition of the Fertil matter
	 \param BurnUp : Maximum achievable burn up envisaged
	 */
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);
	//}
	
	private :
	
	TTree* CreateTMVAInputTree(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);//!<Create input tmva tree to be read by ExecuteTMVA
	double ExecuteTMVA(TTree* theTree);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree


	string fTMVAWeightPath;;//!<The weight needed by TMVA to construct and execute the multilayer perceptron

};

#endif

