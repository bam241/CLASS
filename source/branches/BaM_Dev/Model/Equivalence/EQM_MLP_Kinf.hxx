#ifndef _EQM_MLP_Kinf_HXX
#define _EQM_MLP_Kinf_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"
#include <map>

/*!
 \file
 \brief Header file for EQM_MLP_Kinf class.


 @author BLG
 @version 1.0
 */


using namespace std;

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on neural network to predict @f$k_{\infty}@f$
/*!
The aim of these class is to constuct a fuel from an equivalence model
based on a  Multi layer perceptron (MLP). 
This MLP aims to predict the @f$k_{\infty}(t)@f$ of a PWR-MOX from a given fresh fuel 
composition.
With this MLP prediction and a given number of batch (for the loading plan) an 
average @f$<k_{\infty}>(t)@f$ is calculated according :

@f$<k_{\infty}>^{batch}(t) = \frac{1}{N}\sum_{i}^{N}k_{\infty}(t+\frac{iT}{N})@f$
The maximal reachable burnup has to verify the following conditions :
@f$<k_{\infty}>^{batch}(T/N) = <k_{\infty}>^{batch}(2T/N)  = ... = k_{Threshold}@f$
Where @f$k_{Threshold}@f$ is the criticality threshold which take into account leakage and capture
in non simulated devices such as control rods and mixing grid.

 @author BLG
 @version 1.0
 */
//________________________________________________________________________

class EQM_MLP_Kinf;
#ifndef __CINT__
typedef void (EQM_MLP_Kinf::*PWR_MLP_KINF_DMthPtr)( const string & ) ;
#endif


class EQM_MLP_Kinf : public EquivalenceModel
{
	public :
	/*!
	 \name Constructor
	 */
	//@{

	//{
	/// Polynnomial 2nd order constructor @f$k_{\infty} = \alpha_{0} + \alpha_{1}t + \alpha_{2}t^{2}@f$
 	/// @f$\alpha_{0}@f$, @f$\alpha_{1}@f$, @f$\alpha_{2}@f$ are predict by 3 MLP (one for each)
 	/*!
	 Create a EQM_MLP_Kinf 
	 \param  TMVAWeightPath0 :  PAth to the .xml file containing neural network informations for @f$\alpha_{0}@f$ prediction : PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  TMVAWeightPath1 :  PAth to the .xml file containing neural network informations for @f$\alpha_{1}@f$ prediction: PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  TMVAWeightPath2 :  PAth to the .xml file containing neural network informations for @f$\alpha_{2}@f$ prediction: PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo
	 \param  NumOfBatch : Number of batch for the loading plan (often 3 or 4 for PWR)
	 \param  CriticalityThreshold : Threshold for the @f$k_{\infty}@f$ (see detailed description)
	 */
	EQM_MLP_Kinf(string WeightPathAlpha0, string WeightPathAlpha1, string WeightPathAlpha2, string InformationFile,  int NumOfBatch, double CriticalityThreshold = 1.01);
	//}
	//{
	/// Polynnomial 2nd order constructor @f$k_{\infty} = \alpha_{0} + \alpha_{1}t + \alpha_{2}t^{2}@f$
 	/// @f$\alpha_{0}@f$, @f$\alpha_{1}@f$, @f$\alpha_{2}@f$ are predict by 3 MLP (one for each)
	 /*!
	 Create a EQM_MLP_Kinf 
	 \param log : use for log
	 \param  TMVAWeightPath0 :  PAth to the .xml file containing neural network informations for @f$\alpha_{0}@f$ prediction : PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  TMVAWeightPath1 :  PAth to the .xml file containing neural network informations for @f$\alpha_{1}@f$ prediction: PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  TMVAWeightPath2 :  PAth to the .xml file containing neural network informations for @f$\alpha_{2}@f$ prediction: PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo
	 \param  NumOfBatch : Number of batch for the loading plan (often 3 or 4 for PWR)
	 \param  CriticalityThreshold : Threshold for the @f$k_{\infty}@f$ (see detailed description)
	 */
	EQM_MLP_Kinf(CLASSLogger* log, string WeightPathAlpha0, string WeightPathAlpha1, string WeightPathAlpha2, string InformationFile,  int NumOfBatch, double CriticalityThreshold = 1.01 );
	//}
	//{
	/// Neural network predictor. The kinf(t) is predicted with a MLP 
	/*!
	 Create a EQM_MLP_Kinf 
	 \param  TMVAWeightPath :  PAth to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  NumOfBatch : Number of batch for the loading plan (often 3 or 4 for PWR)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo
	 \param  CriticalityThreshold : Threshold for the @f$k_{\infty}@f$ (see detailed description)
	 */
	EQM_MLP_Kinf(string TMVAWeightPath,int NumOfBatch, string InformationFile = "", double CriticalityThreshold = 1.01);
	//}
	
	//{
	/// Neural network predictor. The kinf(t) is predicted with a MLP 
	/*!
	 Create a EQM_MLP_Kinf
	 \param log : use for log
	 \param  TMVAWeightPath :  PAth to the .xml file containing neural network informations : PATH/TMVAWeight.xml (total path to tmva weight)
	 \param  NumOfBatch : Number of batch for the loading plan (often 3 or 4 for PWR)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo
	 \param  CriticalityThreshold : Threshold for the @f$k_{\infty}@f$ (see detailed description)
	 */
	EQM_MLP_Kinf(CLASSLogger* log, string TMVAWeightPath,int NumOfBatch, string InformationFile = "", double CriticalityThreshold = 1.01);
	//}
	//@}
	
	//{
	/// Return the molar fissile fraction according fissile & ferile content 
	/*!
	 \param Fissil : The composition of the fissile matter
	 \param Fertil : The composition of the Fertil matter
	 \param BurnUp : Maximum achievable burn up envisaged
	 */
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);
	//}
	
	/*!
	 \name Get/Set methods
	 */
	//@{

	void SetBurnUpPrecision(double prop){fBurnUpPrecision = prop;} //!< Set the precision on Burnup : proportion of the targeted burnup
	void SetPCMprecision(double pcm){fPCMprecision = pcm;}		  //!< Set the precision on <k> prediction [pcm]. Neural network predictor constructors

	double GetBurnUpPrecision(){return fBurnUpPrecision;}//!< Get the precision on Burnup : proportion of the targeted burnup
	double GetPCMprecision(){return fPCMprecision/1e5;}//!< Get the precision on <k> prediction []. Neural network predictor constructors
	double GetMaximumBurnUp_MLP(IsotopicVector TheFuel, double TargetBU);

	//@}
	TTree* CreateTMVAInputTree(IsotopicVector FreshFuel, double ThisTime);//!<Create input tmva tree to be read by ExecuteTMVA
	double ExecuteTMVA(TTree* theTree, string WeightPath, bool IsTimeDependant);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree
	double SecondToBurnup(double Second){return Second*fSpecificPower/(24*3.6e6);}

	
	/*!
	 \name Reading NFO related Method
	 */
	//@{
	
	//{
	/// LoadKeyword() : make the correspondance between keyword and reading method
	void LoadKeyword();
	//}
	
	//{
	/// ReadZAIName : read the zai name in the TMWA MLP model
	/*!
	 \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
	 */
	void ReadZAIName(const string &line);
	//}
	
	//{
	/// ReadMaxBurnUp : read a guessed (very overestimated) maximum burnup a fuel can reach (purpose : algorithm initialization)
	/*!
	 \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
	 */
	void ReadMaxBurnUp(const string &line);
	//}
	
	//{
	/// ReadMaxFisContent  : read a guessed (very overestimated) maximum fissile content (purpose : algorithm initialization)
	/*!
	 \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
	 */
	void ReadMaxFisContent(const string &line);
	//}
	
	//{
	/// ReadLine : read a line
	/*!
	 \param line : line to read
	 */
	void ReadLine(string line);
	//}
	
	//@}
	
	private :
	vector <string> fTMVAWeightPath;		//!<The weight needed by TMVA to construct and execute the multilayer perceptron

#ifndef __CINT__
	map<string, PWR_MLP_KINF_DMthPtr> fDKeyword;
#endif
	
	map<ZAI,string> fMapOfTMVAVariableNames;//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step

	
	
	double GetMaximumBurnUp_Pol2(IsotopicVector TheFuel, double TargetBU);
	double GetMaximumBurnUp(IsotopicVector TheFuel, double TargetBU);//!<Get the maximum reachable burnup according the freshfuel composition
	void   GetModelInformation();//!<Read the fMLPInformationFile and fill containers and variables
	
	double BurnupToSecond(double BurnUp){return BurnUp/fSpecificPower*(24*3.6e6);}

	int 	fNumberOfBatch;		//!< The number of batches for the loading plan
	double 	fKThreshold;		//!< The @f$k_{Threshold}@f$
	double 	fSpecificPower; 	//!< The specific power in W/gHM (HM: heavy Metal)
	double 	fMaximalBU;			//!< The approx. maximum burnup reachable by the MLP model		
	double  fMaximalContent;	//!< The approx. maximum fissile content reachable by the MLP model
	double 	fBurnUpPrecision;	//!< precision on Burnup 
	double 	fPCMprecision;		//!< precision on <k> prediction [pcm]
			

};

#endif

