#ifndef _EQUIVALENCEMODEL_
#define _EQUIVALENCEMODEL_


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
#include "CLASSObject.hxx"


using namespace std;

class EquivalenceModel;
#ifndef __CINT__
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
	EquivalenceModel();			//!< Default constructor
	EquivalenceModel(CLASSLogger* log);	//!< Logger constructor
	
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
	virtual	 map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListMassFractionMin, map < string , double> StreamListMassFractionMax, map < string , int> StreamListPriority, map < string , bool> StreamListIsBuffer);
	//}
		
	//@}
	
	/*!
	 \name Get/Set Method
	 */
	//@{

	IsotopicVector GetStreamList(string keyword) {return fStreamList[keyword];}	//!<return the list of ZAI of stream type keyword
	map < string, IsotopicVector> GetAllStreamList() {return fStreamList;}	 	//!<return all the lists

	int GetStreamListNumber(){return fStreamList.size();};
	
	map < string, double> GetBuildFuelFirstGuess(){return fFirstGuessContent;} 	//!<Get the initialization value for BuildFuel algorithm
	map < string, double> GetActualMassContent() const  { return  fActualMassContentInFuel ; }
	map < string, double> GetActualMolarContent() const { return fActualMolarContentInFuel ; }

	virtual map < string , double> GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp) = 0; //!< Return the molar fractions of each element in the fuel according to the burnup, and a given fuel composition (this is the heart of the equivalence model)

	double GetRelativMassPrecision() const	{ return fRelativMassPrecision; }	//!< Mass precision
	int GetMaxInterration()		 const	{ return fMaxInterration; }		//!< Max iterration in build fueld algorythm
	
	double GetBurnUpPrecision(){return fBurnUpPrecision;}//!< Get the precision on Burnup : proportion of the targeted burnup
	
	void SetBuildFuelFirstGuess(map < string, double> FirstGuess){fFirstGuessContent = FirstGuess;} //!<set the initialization value for BuildFuel algorithm (one FistGuess for each component of the fresh fuel
	void SetRelativMassPrecision( double val)	{ fRelativMassPrecision = val; }	//!< Mass precision
	void SetMaxInterration(int val)			{ fMaxInterration = val; }	//!< Max iteration in build fuel algorithm

	
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

	void PrintInfo(); //Print the information red in the INFO stream	
	
	//{
	/// ReadFissil : read the zai name in the TMWA MLP model starts with "k_fissil" keyword
	/*!
	 \param line : line suppossed to contain the fissil list
	 */
	void ReadList(const string &line);

	void ReadFirstGuessContent(const string &line);	
	
	bool isIVInDomain(IsotopicVector IV);

	
	
	
	protected :
	
	map < string, IsotopicVector> fStreamList;					//!< contains all lists of zai needed to build a fuel (example : 2 -> fissileList+fertileList)
											//!< each list is identified by a keyword (example : -> "Fissile" & "Fertile")
	map < string , double>	fStreamListEqMMassFractionMax;		//!< Map that contains lists of stream according to the EqModel with mass maximum fraction
	map < string , double>	fStreamListEqMMassFractionMin;		//!< Map that contains lists of stream according to the EqModel with mass minimum fraction
								
	map < string, double> fFirstGuessContent;					//!< fissile content for BuildFuel initialization (in weight proportion)
	map < string, double> fActualMolarContentInFuel; 				//!< Molar Content in fuel of each list at this step of the calculation
	map < string, double> fActualMassContentInFuel; 				//!< Mass  Content in fuel of each list at this step of the calculation
	double 	fSpecificPower; 							//!< The specific power in W/gHM (HM: heavy Metal)
	double  fMaximalContent;							//!< The approx. maximum fissile content reachable by the model

	
	double LambdaCalculation(string MaterialDenomination, double LambdaPreviousStep, double MaterialMassNeeded, double DeltaMass, vector <IsotopicVector>  Stocks); //!< Calculate the proportion of each stocks in StockArray to take in oder to get a mass of MassNeeded (can be Fer(fertile) or Fis(Fissile))
	void SetLambda(vector<double>& lambda , double Lambda_tot);	//!< Set individual lambda according to the LAMBDA_TOT (lambda of all stocks)	
	void SetLambdaToErrorCode(vector<double>& lambda);
	void StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks);
	double fRelativMassPrecision;		//!< Mass precision
	double 	fBurnUpPrecision;		//!< precision on Burnup 
	int fMaxIterration;			//!< Max iterration in build fueld algorythm

	
	//@}

#ifndef __CINT__
	map<string, EQM_MthPtr> fKeyword;
#endif
	
	bool freaded;
	map< ZAI, pair<double,double> > fZAILimits; 	//!< Fresh fuel range : map<ZAI<min edge ,max edge >>

	string fInformationFile;					//!<  file containing Reactor Type, Fuel type, HM mass, Power, time vector, and TMVA input variables names (looks the manual for format details)
	string fDBFType;					//!<  Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
	string fDBRType;					//!<  Reactor Type (e.g PWR, FBR-Na, ADS..)
	
	/*!
	 \name Others 
	 */
	//@{
	map <string , double > fTotalMassInStocks;  		//!< Total mass in each vector of stock
	map <string , double > fLambdaMax;  		//!< Total lambda of available stocks




	//@}
};

#endif









