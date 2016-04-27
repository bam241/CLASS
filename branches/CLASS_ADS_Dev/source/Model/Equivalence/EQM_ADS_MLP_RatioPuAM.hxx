#ifndef _EQM_ADS_MLP_RatioPuAM_HXX
#define _EQM_ADS_MLP_RatioPuAM_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"

/*!
 \file
 \brief Header file for EQM_ADS_MLP_RatioPuAM

/*!
The aim of these class is to constuct a fuel from a Multi layer perceptron (MLP).
This MLP aims to predict the M.A. fraction such as R = M.A. / (M.A. + Pu) such as it has to verify
 @f$ k(BOC) = k_{target}@f$  with @f$k_{target}@f$ close to 0.96 here.
 WantedTime can't be set by user since it is contain in the .xml file. Indeed this method suppose
 you have trained your MLP to predict the ratio R at BOC.

 @author NT
 @version 1.0
 */
class EQM_ADS_MLP_RatioPuAM : public EquivalenceModel
{
	public :

	EQM_ADS_MLP_RatioPuAM(string TMVAWeightPath, double KeffAtBOC);
	EQM_ADS_MLP_RatioPuAM(CLASSLogger* log, string TMVAWeightPath, double KeffAtBOC);

	virtual map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray);

	void SetSpecificPower(double SpecificPower) {fSpecificPower = SpecificPower;}

	TTree* CreateTMVAInputTree(IsotopicVector TheFuel);
	double ExecuteTMVA(TTree* TheTree, string WeightPath);

	map < string , double> GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp);

	private:

	double BurnupToSecond(double BurnUp){return BurnUp/fSpecificPower*(24*3.6e6);}
	double SecondToBurnup(double Second){return Second*fSpecificPower/(24*3.6e6);}

	int fNumberOfBatch ;		
	double fKThreshold  ;
	double fKeffAtBOC  ;
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

    /**********************************/
    /**	Opening Root file *************/
    /**********************************/
    


};

#endif
