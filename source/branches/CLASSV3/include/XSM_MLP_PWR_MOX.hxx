
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
Design for a PWR MOX reactor

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
	XSM_MLP_PWR_MOX(string TMVA_Weight_Directory,string InformationFile="",bool IsTimeStep=true);
	XSM_MLP_PWR_MOX(CLASSLogger* Log,string TMVA_Weight_Directory,string InformationFile="",bool IsTimeStep=true);	
	~XSM_MLP_PWR_MOX(); 
	//{


 	EvolutionData GetCrossSections(IsotopicVector IV,double t=0) ;


 private :
	void GetDataBaseInformation(); //<! Read information file and fill HM mass, Power, time vector
	vector<double> GetMLPTime() {return fMLP_Time; }

 	void GetMLPWeightFiles();
 	void ReadWeightFile(string Filename, int &Z, int &A, int &I, int &Reaction) ;
 	double ExecuteTMVA(string WeightFile, TTree* InputTree);
 	TTree* CreateTMVAInputTree(IsotopicVector isotopicvector,int TimeStep=0);
 	EvolutionData GetCrossSectionsTime(IsotopicVector IV) ;

 	void ReadWeightFileStep(string Filename, int &Z, int &A, int &I, int &Reaction, int &TimeStep) ;
 	EvolutionData GetCrossSectionsStep(IsotopicVector IV) ;



 	vector<double> 	fMLP_Time;
 	vector<string> 	fWeightFiles;
 	string fTMVAWeightFolder;
 	string fMLPInformationFile;
 	double fDataBasePower;
 	double fDataBaseHMMass;
 	string fDataBaseFType;	//<! Reactor Type (e.g PWR, RNR, ADS..)
 	string fDataBaseRType;	//<! Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
 	bool fIsStepTime;

 	map<ZAI,string> fMapOfTMVAVariableNames;


};

#endif

