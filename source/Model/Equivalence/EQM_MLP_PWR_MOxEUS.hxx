#ifndef _EQM_MLP_PWR_MOxEUS_HXX
#define _EQM_MLP_PWR_MOxEUS_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"

/*!
 \file
 \brief Header file for EQM_MLP_PWR_MOxEUS


 @author FC
 @version 1.0
 */
class EQM_MLP_PWR_MOxEUS : public EquivalenceModel
{
	public :

	EQM_MLP_PWR_MOxEUS(string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold);
	EQM_MLP_PWR_MOxEUS(CLASSLogger* log, string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold);

	void SetSpecificPower(double SpecificPower) {fSpecificPower = SpecificPower;}
	void SetPCMPrecision(double prop) {fPCMPrecision= prop;}	

	double GetPCMPrecision(){return fPCMPrecision/1e5;}

	TTree* CreateTMVAInputTree(IsotopicVector TheFuel, double ThisTime);
	double ExecuteTMVA(TTree* theTree, string WeightPath);
	double CalculateTargetParameter(IsotopicVector FuelToTest);

	private:

	double BurnupToSecond(double BurnUp){return BurnUp/fSpecificPower*(24*3.6e6);}
	double SecondToBurnup(double Second){return Second*fSpecificPower/(24*3.6e6);}

	int fNumberOfBatch ;		
	double fKThreshold  ;
	double fSpecificPower ;	
	double fMaximalBU;
	double 	fPCMPrecision;		

	vector <string> fTMVAWeightPath;		//!<The weight needed by TMVA to construct and execute the multilayer perceptron

	


};

#endif