#ifndef _EQM_FBR_MLP_Keff_HXX
#define _EQM_FBR_MLP_Keff_HXX

#include "EquivalenceModel.hxx"
#include "TTree.h"
#include "TGraph.h"
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

using namespace std;

class EQM_FBR_MLP_Keff;
#ifndef __CINT__
typedef void (EQM_FBR_MLP_Keff::*FBR_MLP_Keff_DMthPtr)( const string & ) ;
#endif

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on neural network to predict @f$k_{eff}@f$.

/*!
The aim of these class is to constuct a fuel from an equivalence model
based on a  Multi layer perceptron (MLP).
This MLP aims to predict the Pu content such as it has to verify
 @f$ k(WantedTime) = k_{target}@f$  with @f$k_{target}@f$ is often close to 1.0 
 but can be set by user. The wanted time is often either 
 the begining of cycle or end of cycle. WantedTime can't be set by user since it is
 contain in the .xml file. Indeed this method suppose you have trained your MLP to predict 
 the @f$k_{eff}@f$ either at BOC or EOC (or any other time)

 @author BLG
 @author BaM
 @version 1.0
 */
//________________________________________________________________________

class EQM_FBR_MLP_Keff : public EquivalenceModel
{
	public:
	/*!
	 \name Constructor
	 */
	//@{



	//{
	///  Create a EQM_FBR_MLP_Keff using Keffective at a given time
 	/// (see class desctiption)
	 /*!
	 Create a EQM_FBR_MLP_Keff 
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations for prediction of keff(t = tfixed)
	 \param  keff_target : Wanted  @f$k_{eff}@f$ (see detailed description)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo (see manual for format)
	 */
	EQM_FBR_MLP_Keff(string TMVAWeightPath, double keff_target = 1.00, string InformationFile = "");
	//}

	//{
	///  Create a EQM_FBR_MLP_Keff using Keffective at a given time
 	/// (see class desctiption)
	 /*!
	 Create a EQM_FBR_MLP_Keff 
	 \param  log:  CLASSLogger object to handle log messages
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations for prediction of keff(t = tfixed)
	 \param  keff_target : Wanted  @f$k_{eff}@f$ (see detailed description)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo (see manual for format)
	 */
	EQM_FBR_MLP_Keff(CLASSLogger* log, string TMVAWeightPath, double keff_target = 1.00, string InformationFile = "");
	//}

	//@}
	
	//{
	/// Return the molar fissile fraction according fissile & ferile content using @f$<k_{\infty}>(t)@f$ prediction
	/*!
	 \param Fissil : The composition of the fissile matter
	 \param Fertil : The composition of the Fertil matter
	 \param BurnUp : Maximum achievable burn up envisaged
	 */
	virtual double GetFissileMolarFraction(vector <IsotopicVector> IVStream, double BurnUp = 0);
	//}
	
	/*!
	 \name Get/Set methods
	 */
	//@{
	void 	SetPCMprecision(double pcm){fPCMprecision = pcm;}	//!< Set the precision on @f$\langle k \rangle@f$ prediction [pcm]. Neural network predictor constructors
	double 	GetPCMprecision(){return fPCMprecision/1e5;}		//!< Get the precision on @f$\langle k \rangle@f$ prediction []. Neural network predictor constructors

	//@}
	
	
	
	
	
	/*!
	 \name TMVA related  methods
	 */
	//@{
	TTree* 	CreateTMVAInputTree(IsotopicVector FreshFuel, double ThisTime);//!<Create input tmva tree to be read by ExecuteTMVA
	
	double 	ExecuteTMVA(TTree* theTree, bool IsTimeDependant);//!<Execute the MLP according to the input tree created by CreateTMVAInputTree
	//@}
	
	
	
	/*!
	 \name Reading NFO related Method
	 */
	//@{
	
	//{
	/// LoadKeyword() : make the correspondance between keyword and reading method
	void LoadKeyword();
	//}
	
	//{
	/// ReadTimeSteps : read the time step of the model
	/*!
	 \param line : line suppossed to contain the time step information starts with "k_timestep" keyword
	 */
	void ReadTimeSteps(const string &line);
	//}
	
	
	//{
	/// ReadZAIName : read the zai name in the TMWA MLP model
	/*!
	 \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
	 */
	void ReadZAIName(const string &line);
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

	string fTMVAWeightPath;					//!<The weight needed by TMVA to construct and execute the multilayer perceptron

#ifndef __CINT__
	map<string, FBR_MLP_Keff_DMthPtr> fDKeyword;
#endif


	map<ZAI,string> fMapOfTMVAVariableNames;//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step

	vector<double> fMLP_Time;	//!< Time (in seconds) when the MLP(t) = keff(t) has been trained.
	double 	fSpecificPower; 	//!< The specific power in W/gHM (HM: heavy Metal)
	double  fMaximalContent;	//!< The approx. maximum fissile content reachable by the MLP model
	

	
	int 	fNumberOfBatch;		//!< The number of batches for the loading plan
	
	double 	fKThreshold;		//!< The @f$k_{Threshold}@f$
	double 	fPCMprecision;		//!< precision on @f$\langle k \rangle@f$ prediction [pcm]
	
	double 	fTargetKeff;		//!< Use for Varying Fissile content to reach fTargetKeff at time used in the MLP Training
	

	
	double 	GetKeffAtFixedTime(IsotopicVector FreshFuel){return ExecuteTMVA( CreateTMVAInputTree(FreshFuel,-1), false );} //!<time independant since the MLP is trained for 1 time

	TGraph* BuildKeffGraph(IsotopicVector FreshFuel);
	TGraph* BuildAverageKeffGraph(TGraph* GRAPH_KEFF);

	double 	GetKeffAt(TGraph* GRAPH_KEFF, int Step);



	
};

#endif

