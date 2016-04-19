#ifndef _EQM_PWR_FixedContent_HXX
#define _EQM_PWR_FixedContent_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"
#include <map>

using namespace std;

class EQM_PWR_FixedContent;
#ifndef __CINT__
typedef void (EQM_PWR_FixedContent::*PWR_Fixed_DMthPtr)( const string & ) ;
#endif

/*!
 \file
 \brief Header file for EQM_PWR_FixedContent


 @author FC
 @version 1.0
 */
class EQM_PWR_FixedContent : public EquivalenceModel
{
	public :

	EQM_PWR_FixedContent(string TMVAWeightPath, int NumOfBatch, string InformationFile = "", double CriticalityThreshold=1.01);
	EQM_PWR_FixedContent(CLASSLogger* log, string TMVAWeightPath, int NumOfBatch, string InformationFile = "", double CriticalityThreshold=1.01);

	virtual	 map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray);

	void SetSpecificPower(double SpecificPower) {fSpecificPower = SpecificPower;}
	void SetBurnUpPrecision(double precision) {fBurnUpPrecision= precision;}
	void SetPCMPrecision(double prop) {fPCMPrecision= prop;}	

	map < string, double> GetFixedMassContent(){return fFixedMassContent;}
	double GetBurnUpPrecision(){return fBurnUpPrecision;}
	double GetPCMPrecision(){return fPCMPrecision/1e5;}

	TTree* CreateTMVAInputTree(IsotopicVector TheFuel, double ThisTime);
	double ExecuteTMVA(TTree* theTree, string WeightPath);
	double GetMaximumBurnUp(IsotopicVector TheFuel, double TargetBU);

	void GetModelInformation();//!<Read the fMLPInformationFile and fill containers and variables
	void LoadKeyword();
	void ReadZAIName(const string &line);
	void ReadMaxBurnUp(const string &line);
	void ReadMaxFisContent(const string &line);
	void ReadFixedMassContent(const string &line);
	void ReadLine(string line);

	map < string , double> GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp);

	double BurnupToSecond(double BurnUp){return BurnUp/fSpecificPower*(24*3.6e6);}
	double SecondToBurnup(double Second){return Second*fSpecificPower/(24*3.6e6);}

	private :
	vector <string> fTMVAWeightPath;		//!<The weight needed by TMVA to construct and execute the multilayer perceptron

#ifndef __CINT__
	map<string, PWR_Fixed_DMthPtr> fDKeyword;
#endif

	int fNumberOfBatch ;		
	double fKThreshold  ;
	double fMaximalBU;
	double fBurnUpPrecision;
	double fPCMPrecision;	
	map < string, double> fFixedMassContent;	
	
	map<ZAI,string> fMapOfTMVAVariableNames;//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step	

};

#endif