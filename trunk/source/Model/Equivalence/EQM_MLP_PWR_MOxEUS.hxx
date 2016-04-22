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

	EQM_MLP_PWR_MOxEUS(string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold, double MaximalPuMassContent);
	EQM_MLP_PWR_MOxEUS(CLASSLogger* log, string TMVAWeightPath, int NumOfBatch, double CriticalityThreshold, double MaximalPuMassContent);

	virtual	 map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray);

	void SetSpecificPower(double SpecificPower) {fSpecificPower = SpecificPower;}
	void SetMinimalU5Enrichment(double MinEU5) {fMinimalU5Enrichment= MinEU5;}
	void SetMaximalU5Enrichment(double MaxEU5){fMaximalU5Enrichment= MaxEU5;}
	void SetBurnUpPrecision(double precision) {fBurnUpPrecision= precision;}
	void SetPCMPrecision(double prop) {fPCMPrecision= prop;}	

	double GetBurnUpPrecision(){return fBurnUpPrecision;}
	double GetPCMPrecision(){return fPCMPrecision/1e5;}
	double GetMinimalU5Enrichment(){return fMinimalU5Enrichment;};

	TTree* CreateTMVAInputTree(IsotopicVector TheFuel, double ThisTime);
	double ExecuteTMVA(TTree* theTree, string WeightPath);
	double GetMaximumBurnUp(IsotopicVector TheFuel, double TargetBU);
	double GetU5Enrichment(map <string , IsotopicVector > IVStream, double TargetBU);

	map < string , double> GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp);

	private:

	double BurnupToSecond(double BurnUp){return BurnUp/fSpecificPower*(24*3.6e6);}
	double SecondToBurnup(double Second){return Second*fSpecificPower/(24*3.6e6);}

	int fNumberOfBatch ;		
	double fKThreshold  ;
	double fSpecificPower ;	
	double fMaximalBU;
	double fBurnUpPrecision;
	double 	fPCMPrecision;		
	
	double fActualPuMolarContent;
	double fActualPuMassContent;		
	double fMaximalPuMassContent ;
	double fMaximalPuMolarContent;
	double fMinimalU5Enrichment;
	double fMaximalU5Enrichment;

	vector <string> fTMVAWeightPath;		//!<The weight needed by TMVA to construct and execute the multilayer perceptron

	


};

#endif