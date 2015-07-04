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

/*!
 \file
 \brief Header file for EQM_FBR_MLP_Keff class.


 @author BLG
 @version 1.0
 */


using namespace std;

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on neural network to predict @f$k_{eff}@f$
/*!
The aim of these class is to constuct a fuel from an equivalence model
based on a  Multi layer perceptron (MLP). 
This MLP aims to predict either :

The Pu content is set such as it has to verify
 @f$ k(WantedTime) = k_{target}@f$  with @f$k_{target}@f$ is often close to 1.0 
 but can be set by user. The wanted time is often either 
 the begining of cycle or end of cycle. WantedTime can't be set by user since it is
 contain in the .xml file. Indeed this method suppose you have trained your MLP to predict 
 the keffective either at BOC or EOC (or any other time)

 \warning 
 With method 1, it is not guaranted that there is a solution for Pu content verifying : 
 @f$<k_{\infty}>^{batch}(t) = \frac{1}{N}\sum_{i}^{N}k_{\infty}(t+\frac{iT}{N} )@f$

 @author BLG
 @author BaM
 @version 1.0
 */
//________________________________________________________________________
class EQM_FBR_MLP_Keff;
#ifndef __CINT__
typedef void (EQM_FBR_MLP_Keff::*MLP_FBR_Keff_DMthPtr)( const string & ) ;
#endif


class EQM_FBR_MLP_Keff : public EquivalenceModel
{
	public :
	/*!
	 \name Constructor
	 */
	//@{



	//{
	///  Create a EQM_FBR_MLP_Keff using Keffective at a given time
 	/// (see class desctiption)
	 /*!
	 Create a EQM_FBR_MLP_Keff 
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations for prediction of keff(t=tfixed)
	 \param  keff_target : Wanted  @f$k_{eff}@f$ (see detailed description)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo (see manual for format)
	 */
	EQM_FBR_MLP_Keff(string TMVAWeightPath, double keff_target = 1.00, string InformationFile="");
	//}

	//{
	///  Create a EQM_FBR_MLP_Keff using Keffective at a given time
 	/// (see class desctiption)
	 /*!
	 Create a EQM_FBR_MLP_Keff 
	 \param  log:  CLASSLogger object to handle log messages
	 \param  TMVAWeightPath :  Path to the .xml file containing neural network informations for prediction of keff(t=tfixed)
	 \param  keff_target : Wanted  @f$k_{eff}@f$ (see detailed description)
	 \param  InformationFile : Total path to the file containing time steps, fissile and ferile list (ante and post fabrication time cooling). Default is the same total path as TMVAWeightPath but extension is replaced by .nfo (see manual for format)
	 */
	EQM_FBR_MLP_Keff(CLASSLogger* log, string TMVAWeightPath, double keff_target = 1.00, string InformationFile="");
	//}

	//@}
	
	//{
	/// Return the molar fissile fraction according fissile & ferile content using @f$<k_{\infty}>(t)@f$ prediction
	/*!
	 \param Fissil : The composition of the fissile matter
	 \param Fertil : The composition of the Fertil matter
	 \param BurnUp : Maximum achievable burn up envisaged
	 */
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp=0);
	//}
	
	/*!
	 \name Get/Set methods
	 */
	//@{
	void 	SetPCMprecision(double pcm){fPCMprecision = pcm;}	//!< Set the precision on <k> prediction [pcm]. Neural network predictor constructors
	double 	GetPCMprecision(){return fPCMprecision/1e5;}		//!< Get the precision on <k> prediction []. Neural network predictor constructors

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
	/// ReadSpecificPower : read the Specific Power of the DataBase
	/*!
	 \param line : line suppossed to contain the Specific Power information starts with "k_specpower" keyword
	 */
	void ReadSpecificPower(const string &line);
	//}

	//{
	/// ReadMaximalContent : read the approx. maximum fissile content reachable by the MLP model
	/*!
	 \param line : line suppossed to contain the maximal content information starts with "k_contentmax" keyword
	 */
	void ReadMaximalContent(const string &line);
	//}
	
	//{
	/// ReadZAIName : read the zai name in the TMWA MLP model
	/*!
	 \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
	 */
	void ReadZAIName(const string &line);
	//}
	
	//{
	/// ReadFissil : read the zai name in the TMWA MLP model starts with "k_fissil" keyword
	/*!
	 \param line : line suppossed to contain the fissil list
	 */
	void ReadFissil(const string &line);
	//}

	//{
	/// ReadFertil : read the zai name in the TMWA MLP model starts with "k_fertil" keyword
	/*!
	 \param line : line suppossed to contain the fertil list & their default isotopic fraction
	 */
	void ReadFertil(const string &line);
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
	map<string, MLP_FBR_Keff_DMthPtr> fDKeyword;
#endif

	void   	GetModelInformation();//!<Read the fMLPInformationFile and fill containers and variables
	
	map<ZAI,string> fMapOfTMVAVariableNames;//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step

	vector<double> fMLP_Time;	//!< Time (in seconds) when the MLP(t)=keff(t) has been trained.
	double 	fSpecificPower; 	//!< The specific power in W/gHM (HM: heavy Metal)
	double  fMaximalContent;	//!< The approx. maximum fissile content reachable by the MLP model
	

	
	
	
	int 	fNumberOfBatch;		//!< The number of batches for the loading plan
	
	double 	fKThreshold;		//!< The @f$k_{Threshold}@f$
	double 	fPCMprecision;		//!< precision on <k> prediction [pcm]
	
	bool	fIsAverageKeff;		//!< True if using the first contructor (fissile content prediction with average keff (over batches))
	double  fKmin;				//!< Lower edge of kedd Used by second constructor (fissile content prediction using keff at BOC (or other time)
	double  fKmax;				//!< Upper edge of kedd Used by second constructor (fissile content prediction using keff at BOC (or other time)
	double 	fTargetKeff;		//!< Use for Varying Fissile content to reach fTargetKeff at time used in the MLP Training
	

	
	

	
	double 	GetKeffAtFixedTime(IsotopicVector FreshFuel){return ExecuteTMVA( CreateTMVAInputTree(FreshFuel,-1), false );} //!<time independant since the MLP is trained for 1 time

	TGraph* BuildKeffGraph(IsotopicVector FreshFuel);
	TGraph* BuildAverageKeffGraph(TGraph* GRAPH_KEFF);

	double 	GetKeffAt(TGraph* GRAPH_KEFF, int Step);



	
};

#endif

