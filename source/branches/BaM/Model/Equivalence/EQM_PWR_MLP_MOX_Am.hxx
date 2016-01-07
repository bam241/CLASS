#ifndef _EQM_PWR_MLP_MOX_AM_HXX
#define _EQM_PWR_MLP_MOX_AM_HXX

#include "EquivalenceModel.hxx"
#include "TMVA/Reader.h"

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
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 */
	EQM_PWR_MLP_MOX_AM(string TMVAWeightPath);
	//}
	
	//{
	/// Logger constructor
	/*!
	 Create a EQM_PWR_MLP_MOX_AM
	 \param log : use for log
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
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
	
	/*!
	 \name TMVA related methods
	 */
	//@{

		void UpdateInputComposition(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);//!<Create input tmva tree to be read by ExecuteTMVA
		double ExecuteTMVA(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree

	//@}	

	private :

    float Pu8;
    float Pu9;
    float Pu10;
    float Pu11;
    float Pu12;
    float Am1;
    float Am2;
    float Am3;

    float BU;
    float U5_enrichment;

    TMVA::Reader *reader;
	string fTMVAWeightPath;;//!<The weight needed by TMVA to construct and execute the multilayer perceptron

};

#endif

