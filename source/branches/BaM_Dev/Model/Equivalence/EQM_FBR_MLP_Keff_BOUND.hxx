#ifndef _EQM_FBR_MLP_Keff_BOUND_HXX_
#define _EQM_FBR_MLP_Keff_BOUND_HXX_

#include "EquivalenceModel.hxx"
#include "TTree.h"
#include "TGraph.h"

using namespace std;

class EQM_FBR_MLP_Keff_BOUND;
#ifndef __CINT__
typedef void (EQM_FBR_MLP_Keff_BOUND::*FBR_MLP_Keff_BOUND_DMthPtr)( const string & ) ;
#endif

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on neural network to predict @f$k_{\infty}@f$

/*!
 The aim of these class is to constuct a fuel from an equivalence model
 based on a  Multi layer perceptron (MLP).
 This MLP aims to predict :
 
 The @f$k_{\infty}(t)@f$ of a FBR-MOX from a given fresh fuel composition.
 With this MLP prediction and a given number of batch (for the loading plan) an
 average @f$<k_{\infty}(t)>@f$ is calculated according :
 @f$<k_{\infty}>^{batch}(t) = \frac{1}{N}\sum_{i}^{N} k_{\infty}(t + \frac{iT}{N} )@f$
 The Fissile content has to verify this condition :
 @f$ k_{\infty Max} \geq <k_{\infty}>^{batch}(T/N) \geq k_{\infty Min} @f$
 Where @f$ k_{\infty Max}@f$ and @f$k_{\infty Min}@f$ are arguments of the constructor.
 
 
 \warning
 it is not guaranted that there is a solution for Pu content verifying :
 @f$ k_{\infty Max} \geq <k_{\infty}>^{batch}(T/N) \geq k_{\infty Min} @f$
 
 @author BLG
 @author BaM
 @version 1.0
 */
//________________________________________________________________________

class EQM_FBR_MLP_Keff_BOUND : public EquivalenceModel
{
	public:
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	///  Create a EQM_FBR_MLP_Keff_BOUND using @f$k_{\infty}@f$ average over batch
	/// (see class desctiption)
	/*!
	 Create a EQM_FBR_MLP_Keff_BOUND using @f$k_{\infty}@f$ average over batch
	 \param  TMVAWeightPath0:  Path to the .xml file containing neural network informations for prediction of keff(t)
	 \param  NumOfBatch : Number of batch for the loading plan (often 5 for FBR-Na)
	 \param  KeffMin : Lower average @f$<k_{\infty}>@f$ value
	 \param  KeffMax : Upper average @f$<k_{\infty}>@f$ value
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo (see manual for format)
	 */
	EQM_FBR_MLP_Keff_BOUND(string TMVAWeightPath,  int NumOfBatch , double LowerKeffective, double UpperKeffective, string InformationFile = "");
	//}
	
	//{
	///  Create a EQM_FBR_MLP_Keff_BOUND using @f$k_{\infty}@f$ average over batch
	/// (see class desctiption)
	/*!
	 Create a EQM_FBR_MLP_Keff_BOUND using @f$k_{\infty}@f$ average over batch
	 \param  log:  CLASSLogger object to handle log messages
	 \param  TMVAWeightPath0:  Path to the .xml file containing neural network informations for prediction of keff(t)
	 \param  NumOfBatch : Number of batch for the loading plan (often 5 for FBR-Na)
	 \param  KeffMin : Lower average @f$<k_{\infty}>@f$ value
	 \param  KeffMax : Upper average @f$<k_{\infty}>@f$ value
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo (see manual for format)
	 */
	EQM_FBR_MLP_Keff_BOUND(CLASSLogger* log, string TMVAWeightPath,  int NumOfBatch , double LowerKeffective, double UpperKeffective, string InformationFile = "");
	//}
	
	
	//@}
	
	//{
	/// Return the molar fissile fraction according fissile & ferile content using @f$<k_{\infty}>(t)@f$ prediction
	/*!
	 \param Fissil : The composition of the fissile matter
	 \param Fertil : The composition of the Fertil matter
	 \param BurnUp : Maximum achievable burn up envisaged
	 */
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp = 0);
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
	/// ReadLine : read a line
	/*!
	 \param line : line to read
	 */
	void ReadLine(string line);
	//}
	
	//@}
	
	
	
	private :
	
	string fTMVAWeightPath;					//!<The weight needed by TMVA to construct and execute the multilayer perceptron
	string fMLPInformationFile;				//!<The path to the informations necessary to execute the MLP
	
#ifndef __CINT__
	map<string, FBR_MLP_Keff_BOUND_DMthPtr> fDKeyword;
#endif
	
	map<ZAI,string> fMapOfTMVAVariableNames;//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step
	
	vector<double> fMLP_Time;	//!< Time (in seconds) when the MLP(t) = keff(t) has been trained.
	
	
	
	
	
	int 	fNumberOfBatch;		//!< The number of batches for the loading plan
	
	double 	fKThreshold;		//!< The @f$k_{Threshold}@f$
	double 	fPCMprecision;		//!< precision on @f$\langle k \rangle@f$ prediction [pcm]
	double  fKmin;				//!< Lower edge of kedd Used by second constructor (fissile content prediction using keff at BOC (or other time)
	double  fKmax;				//!< Upper edge of kedd Used by second constructor (fissile content prediction using keff at BOC (or other time)
	double 	fTargetKeff;		//!< Use for Varying Fissile content to reach fTargetKeff at time used in the MLP Training
	
	
	
	
	
	/*!
	 \name kinf prediction methods & keff averaging
	 */
	//@{
	 
	double 	GetKeffAtFixedTime(IsotopicVector FreshFuel){return ExecuteTMVA( CreateTMVAInputTree(FreshFuel,-1), false );} //!<time independant since the MLP is trained for 1 time
	
	TGraph* BuildKeffGraph(IsotopicVector FreshFuel);
	TGraph* BuildAverageKeffGraph(TGraph* GRAPH_KEFF);
	
	double 	GetKeffAt(TGraph* GRAPH_KEFF, int Step);
	
	//@}
	
	
	
};

#endif








