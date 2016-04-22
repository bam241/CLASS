#ifndef _EQM_PWR_MLP_MOX_HXX
#define _EQM_PWR_MLP_MOX_HXX

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
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_PWR_MLP_MOX : public EquivalenceModel
{
	public :
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// normal constructor
	/*!
	 Create a EQM_PWR_MLP_MOX 
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 */
	EQM_PWR_MLP_MOX(string TMVAWeightPath);
	//}
	
	//{
	/// Logger constructor
	/*!
	 Create a EQM_PWR_MLP_MOX
	 \param log : use for log
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 */
	EQM_PWR_MLP_MOX(CLASSLogger* log, string TMVAWeightPath);
	//}
	//@}

    ~EQM_PWR_MLP_MOX();

	
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

    void UpdateInputComposition(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);
    double ExecuteTMVA(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree
		
	//@}	

	private :

    float fPu8;
    float fPu9;
    float fPu10;
    float fPu11;
    float fPu12;
    float fAm1;

    float fBU;
    float fU5_enrichment;

    TMVA::Reader *freader;
	string fTMVAWeightPath;;//!<The weight needed by TMVA to construct and execute the multilayer perceptron

};

#endif

