
#ifndef _XSM_MLP_PWR_MOX_HXX
#define _XSM_MLP_PWR_MOX_HXX


/*!
 \file
 \brief Header file for XSM_MLP_PWR_MOX class.
 
 
 @authors BLG
 @version 1.0
 */
#include "XSModel.hxx"
#include "TTree.h"
#include <string>
#include <fstream>
#include <iostream> 
#include <map>
#include <vector>
typedef long long int cSecond;
using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a XSM_MLP_PWR_MOX.
This is the class to predict cross sections with a set of MultiLayerPerceptrons

 @authors BLG
 @version 1.0
 */
//________________________________________________________________________


class XSM_MLP_PWR_MOX : public XSModel
{
public : 

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	/// CONSTRUCTOR 
	/*!
	 \param string TMVA_Weight_Directory The directory where all the TMVA weight are located
	 \param string InformationFile, Name of the information file located in TMVA_Weight_Directory (defeult: Data_Base_Info.nfo)
	 \param bool IsTimeStep, if true , one TMVA weihgt per step time is requiered otherwise it assumes time is part of the MLP inputs

	 */
	XSM_MLP_PWR_MOX(string TMVA_Weight_Directory,string InformationFile="",bool IsTimeStep=true);
	XSM_MLP_PWR_MOX(CLASSLogger* Log,string TMVA_Weight_Directory,string InformationFile="",bool IsTimeStep=false);	
	~XSM_MLP_PWR_MOX(); 
	//{


 	EvolutionData GetCrossSections(IsotopicVector IV,double t=0) ;//!< Return calculated cross section by the MLP regression


 private :
	void GetDataBaseInformation(); //<! Read information file and fill Reactor Type, Fuel type,  HM mass, Power, time vector, and TMVA input variables names (looks at Data Bases example for format details)
	vector<double> GetMLPTime() {return fMLP_Time; }//<! Return time vector (seconds) defined in the DataBaseInformation file

 	void GetMLPWeightFiles();//<! Find all .xml file in TMVA_Weight_Directory
 	void ReadWeightFile(string Filename, int &Z, int &A, int &I, int &Reaction) ;//<! Select the reaction according to the weight file name (file name is formated !!) read the manual to know it (formated for fIsTimeStep==false)
 	double ExecuteTMVA(string WeightFile, TTree* InputTree);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree
 	TTree* CreateTMVAInputTree(IsotopicVector isotopicvector,int TimeStep=0);//!<Create input tmva tree to be read by ExecuteTMVA
 	EvolutionData GetCrossSectionsTime(IsotopicVector IV) ;//!< Return calculated cross section by the MLP regression when fIsTimeStep==false

 	void ReadWeightFileStep(string Filename, int &Z, int &A, int &I, int &Reaction, int &TimeStep) ;//<! Select the reaction according to the weight file name (file name is formated !!) read the manual to know it (formated for fIsTimeStep==true)
 	EvolutionData GetCrossSectionsStep(IsotopicVector IV) ;//!< Return calculated cross section by the MLP regression when fIsTimeStep==true



 	vector<double> 	fMLP_Time;//<! Time vector of the data base
 	vector<string> 	fWeightFiles;//<! All the weight file contains in fTMVAWeightFolder
 	string fTMVAWeightFolder;//<!  folder containing all the weight file 
 	string fMLPInformationFile;//<! file containing Reactor Type, Fuel type, HM mass, Power, time vector, and TMVA input variables names (looks the manual for format details)
 	double fDataBasePower;//<!Power of the data base (read from fMLPInformationFile )
 	double fDataBaseHMMass;//<!Heavy metal mass of the data base (read from fMLPInformationFile )
 	string fDataBaseFType;	//<! Reactor Type (e.g PWR, RNR, ADS..)
 	string fDataBaseRType;	//<! Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
 	bool fIsStepTime;//<!true if one TMVA weihgt per step time is requiered otherwise it assumes time is part of the MLP inputs 

 	map<ZAI,string> fMapOfTMVAVariableNames;//<! List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step


};

#endif

