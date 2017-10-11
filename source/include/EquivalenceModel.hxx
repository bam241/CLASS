#ifndef _FUELLOADINGMODEL_
#define _FUELLOADINGMODEL_


/*!
 \file
 \brief Header file for EquivalenceModel class.


 @author BLG
 @author BaM
 @author FaC
 @version 3.0
 */

#include "IsotopicVector.hxx"
#include <math.h>
#include "TTree.h"
#include <map>
#include "CLASSObject.hxx"


using namespace std;

class EquivalenceModel;
#ifndef __ROOTCLING__
typedef void (EquivalenceModel::*EQM_MthPtr)( const string & );
#endif

//-----------------------------------------------------------------------------//

//! Determines how to build a fresh fuel
/*!
 Define an EquivalenceModel.
 The aim of these class is to gather all the commum properties of all
 Equivalence Model.

 \warning
 Never instantiate EquivalenceModel in your CLASS input but it's derivated class
 @see EQM_FBR_BakerRoss_MOX.
 @see EQM_PWR_MLP_MOX
 @see EQM_FBR_MLP_Keff.hxx
 @see EQM_PWR_MLP_MOX_Am
 @see EQM_FBR_MLP_Keff_BOUND
 @see EQM_PWR_POL_UO2
 @see EQM_MLP_Kinf.hxx
 @see EQM_PWR_QUAD_MOX
 @see EQM_PWR_LIN_MOX

 @author BLG
 @author BaM
 @author FaC

 @version 3.0
 */
//________________________________________________________________________


class EquivalenceModel : public CLASSObject
{
public :
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	EquivalenceModel();			//!< Default constructor with path
	EquivalenceModel(CLASSLogger* log);	//!< Logger constructor with path

	virtual ~EquivalenceModel();		//!< Destructor
	//@}
	
	/*!
	 \name Fuel Construction Method
	 */
	//@{
	
	//{
	/// BuildFuel function.
	/*!
	 Build the fuel following the equivalance model with the proper requierment in term of mass, burnup....
	 \param double burnup reached by the fuel at the end of irradiation
	 \param double HMMass, Heavy metal mass needed
    	 \param map < string , vector <IsotopicVector> > StreamArray, the string is the stream code (fissile fertile ,...) the IsotopicVector the fraction of each IV to take in the (fissile, fertile,..) stock .
     */
	virtual	 map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListMassFractionMin, map < string , double> StreamListMassFractionMax, map < int , string> StreamListPriority, map < string , bool> StreamListIsBuffer);
	//}


	//@}
	TTree* CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)	; //!<Create input tmva tree to be read by ExecuteTMVA
	double CalculateTargetParameter(IsotopicVector TheFuel, string TargetParameterName); //!<Get a fuel parameter associated to the fuel ---> ex : BurnUpMax, keffBOC, keffEOC, ... 
	double CalculateBurnUpMax(IsotopicVector TheFuel, map<string, double> ModelParameter);//!<Calculate the BU max associated to a fuel composition based on MLP prediction (suitable for PWR)
	double CalculateKeffAtBOC(IsotopicVector TheFuel); //!<Calculate the keff at BOC associated to a fuel composition based on MLP prediction (suitable for SFR)

    	string GetTMVAXMLFilePath() {return fTMVAXMLFilePath;} // Return the path to TMVA XML File path
    	string GetTMVANFOFilePath() {return fTMVANFOFilePath;} // Return the path to TMVA NFO File path
    	void SetTMVAXMLFilePath(string TMVAXMLFilePath) {fTMVAXMLFilePath = TMVAXMLFilePath;} // Set the path to TMVA XML File path
    	void SetTMVANFOFilePath(string TMVANFOFilePath) {fTMVANFOFilePath = TMVANFOFilePath;} // Set the path to TMVA NFO File path

    /*!
	 \name Get/Set Method
	 */
	//@{

	IsotopicVector GetStreamList(string keyword) {return fStreamList[keyword];}	//!<return the list of ZAI of stream type keyword
	map < string, IsotopicVector> GetAllStreamList() {return fStreamList;}	 	//!<return all the lists

	int GetStreamListNumber(){return fStreamList.size();};
	int GetMaxIterration()		 const	{ return fMaxIterration; }		//!< Max iterration in build fueld algorythm	
	double GetTargetParameterStDev(){return fTargetParameterStDev;}//!< Get the precision on fTargetParameterStDev
	double GetStreamListEqMMassFractionMax(string keyword){return fStreamListEqMMassFractionMax[keyword] ;}
	double GetStreamListEqMMassFractionMin(string keyword){return fStreamListEqMMassFractionMin[keyword] ;}
	double GetPCMPrecision(){return fPCMprecision/1e5;}//!< Get the precision on @f$\langle k \rangle@f$ prediction []. Neural network predictor constructors

    	void SetModelParameter(string sMP, double dMP)  { fModelParameter[sMP] = dMP; }   //!< Set Model Parameters precised in NFO file
    	map<string, double> GetModelParameter()  { return fModelParameter; }   //!< Get Model Parameters precised in NFO file

    void SetNonZaiTMVAVariable(string snZP, double dnZP);    //!< Set NonZaiTMVAVariables
	vector< pair<double, string> > GetNonZaiTMVAVariables()  { return fListOfNonZaiTMVAVariables; }   //!< Get NonZaiTMVAVariables

	void SetMaxIterration(int val)	{ fMaxIterration = val; }	//!< Max iteration in build fuel algorithm
	void SetTargetParameterStDev(double TPSD){fTargetParameterStDev = TPSD;} //!< Set the precision on Target Parameter
	void SetStreamListEqMMassFractionMax(string keyword, double value){fStreamListEqMMassFractionMax[keyword] = value;}
	void SetStreamListEqMMassFractionMin(string keyword, double value){fStreamListEqMMassFractionMin[keyword] = value;}

	void SetPCMPrecision(double pcm){fPCMprecision = pcm;}		  //!< Set the precision on @f$\langle k \rangle@f$ prediction [pcm]. Neural network predictor constructors

	/*!
	 \name Time <-> Burnup conversion
	 */
	//@{

	double SecondToBurnup(double Second){return Second*fSpecificPower/(24*3.6e6);}
	double BurnupToSecond(double BurnUp){return BurnUp/fSpecificPower*(24*3.6e6);}

	//@}

	//@}
	
	/*!
	 \name Reading NFO related Method
	 */
	//@{
	void ReadNFO();
	virtual void ReadLine(string line);

	
	virtual void LoadKeyword();
	void ReadZAIlimits(const string &line);
	void ReadType(const string &line);
	
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
	 \param line : line suppossed to contain the ZAI name  starts with "k_maxburnup" keyword
	 */
	void ReadMaxBurnUp(const string &line);
	//}
	
	//{
	/// ReadSpecificPower : read the Specific Power of the DataBase
	/*!
	 \param line : line suppossed to contain the Specific Power information starts with "k_specpower" keyword
	 */
	void ReadSpecificPower(const string &line);
	//}

    map < string, IsotopicVector> GetAllStreamList() {return fStreamList;}		// return all the lists

	//{
	/// ReadNonZaiTMVAVariables : read the NonZai variables for the predictor (ex : Nbatch, Specific power)
	/*!
	 \param line : line suppossed to contain the NonZai variables for TMVA starts with "k_nonZAIforTMVA" keyword
	 */
	void ReadNonZaiTMVAVariables(const string &line);
	//}

	//{
	/// ReadOutput : read the output type of the predictor (ex : kinf)
	/*!
	 \param line : line suppossed to contain the Specific Power information starts with "k_output" keyword
	 */
	void ReadOutput(const string &line);
	//}

	bool isIVInDomain(IsotopicVector IV);
	void StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks);
	void ConvertMassToLambdaVector(string MaterialDenomination, vector<double>& lambda, double MaterialMassNeeded, vector <IsotopicVector>  Stocks);

	string fPredictorType ; 								//!< Type of predictor used in Equivalence Model (ex: MLP)
	string fOutput ; 								//!< Type of output calculated by the predictor
	string fBuffer ; 									//!< Name of material used as buffer in fuel 

	map<string, double> fModelParameter ; 					///< Map of equivalence model parameter 

 	vector< pair<double, string> > fListOfNonZaiTMVAVariables ; 					///!<  List of TMVA input variable names that are not ZAIs
 	map<ZAI,string> fMapOfTMVAVariableNames;				//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step

	double 	fTargetParameterStDev;							//!< Precision on target parameter calculation 

    map < string, IsotopicVector> fStreamList;                  //!< contains all lists of zai needed to build a fuel (example : 2 -> fissileList+fertileList)
                                                                //!< each list is identified by a keyword (example : -> "Fissile" & "Fertile")
	double 	fSpecificPower; 					         		//!< The specific power in W/gHM (HM: heavy Metal)
	void SetLambdaToErrorCode(vector<double>& lambda);


#ifndef __ROOTCLING__
	map<string, EQM_MthPtr> fKeyword;
#endif

    bool freaded;

	map< ZAI, pair<double, double> > fZAILimits; 	//!< Fresh fuel range : map<ZAI<min edge ,max edge >>

	/*!
	 \name Others
	 */
	//@{
	map <string , double > fTotalMassInStocks;      //!< Total mass in each vector of stock
	map <string , double > fLambdaMax;  			//!< Total lambda of available stocks

	//@}

private :

};

#endif









