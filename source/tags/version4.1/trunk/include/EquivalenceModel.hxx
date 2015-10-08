#ifndef _EQUIVALENCEMODEL_
#define _EQUIVALENCEMODEL_


/*!
 \file
 \brief Header file for EquivalenceModel class.
 
 
 @author BaM
 @author BLG
 @version 3.0
 */

#include "IsotopicVector.hxx"
#include <math.h>
#include <map>
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
	 \param vector<double> &lambda, fraction of the stock to take (initialy should be 0)
	 \param vector<IsotopicVector> FissilArray, isotopicvectors to use to get the fissile part of the fuel
	 \param vector<IsotopicVector> FertilArray, isotopicvectors to use to get the fertile part of the fuel (if empty take it from the OutIncome)
	 */
	virtual	 vector<double> BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray );
	//}
		
	//@}
	
	/*!
	 \name Get/Set Method
	 */
	//@{
	
	IsotopicVector GetFertileList() {return fFertileList;}	//!<return the fertile list
	IsotopicVector GetFissileList() {return fFissileList;}	//!<return the fissile list
	
	double GetBuildFuelFirstGuess(){return fFirstGuessFissilContent;} //!<Get the initialization value for BuildFuel algorithm
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp) {return 0;} //!< Return the molar fraction of fissile element in the fuel according to the burnup, and a given fuel composition (this is the heart of the equivalence model)

	double GetRelativMassPrecision() const	{ return fRelativMassPrecision; }	//!< Mass precision
	int GetMaxInterration()		 const	{ return fMaxInterration; }		//!< Max iterration in build fueld algorythm

	double GetActualFissileContent() const { return fActualFissileContent; }	//!< Get the fissile content at the actual algorithm step (usefull for EQM_MLP_Kinf)
	
	void SetFertileList(IsotopicVector IV) {fFertileList = IV;}//!<set the fertile list
	void SetFissileList(IsotopicVector IV) {fFissileList = IV;}//!<set the fissile list
	void SetBuildFuelFirstGuess(double FirstGuess){fFirstGuessFissilContent = FirstGuess;} //!<set the initialization value for BuildFuel algorithm
	void SetRelativMassPrecision( double val)	{ fRelativMassPrecision = val; }	//!< Mass precision
	void SetMaxInterration(int val)			{ fMaxInterration = val; }		//!< Max iterration in build fueld algorythm

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
	
	void PrintInfo(); //Print the information red in the INFO stream

	//@}
	
	
	bool isIVInDomain(IsotopicVector IV);

	
	protected :
	
	IsotopicVector fFertileList;	//!< contain the list of zai, needed as fertile, taken in a stock before fabrication
	//!< if no stock are provided, will take the isotopic vector in the Park income
	
	IsotopicVector fFissileList;	//!< contain the list of zai, needed as fissile, taken in a stock before fabrication
	//!< if no stock are provided the fuel will not be made
	
	double fFirstGuessFissilContent;//!< fissile content for BuildFuel initialization (in weight proportion)
	double fActualFissileContent;	//!< fissile content at the actual dichotomy step (usefull for EQM_MLP_Kinf)
	
	
	protected :
	/*!
	 \name Algorithm parameters
	 */
	//@{	
	double 	fRelativMassPrecision;		//!< Mass precision
	int 	fMaxInterration;			//!< Max iterration in build fueld algorythm
	//@}

	/*!
	 \name Algorithm related functions for default BuildFuel function
	 */
	//@{
	void 	SetLambda(vector<double>& lambda ,int FirstStockID, int LastStockID, double LAMBDA_TOT);	//!< Set individual lambda according to the LAMBDA_TOT (lambda of all stocks)
	void	SetLambdaToErrorCode(vector<double>& lambda);												//!< Set all vector elements to -1
	double 	LAMBDA_TOT_FOR(double MassNeeded, vector<IsotopicVector> StockArray, string FisOrFer);		//!< Calculate the proportion of each stocks in StockArray to take in oder to get a mass of MassNeeded (can be Fer(fertile) or Fis(Fissile))
	bool 	Build_Fuel_According_Lambda(vector<double> &lambda,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray, double HMMass,IsotopicVector &Fissile, IsotopicVector &Fertile);	//!< Compute the fuel isotopy according to the proportion of each IV taken
	
	//@}

#ifndef __CINT__
	map<string, EQM_MthPtr> fKeyword;
#endif
	
	bool freaded;
	map< ZAI, pair<double,double> > fZAILimits; //!< Fresh fuel range : map<ZAI<min edge ,max edge >>

	string fInformationFile;	//!<  file containing Reactor Type, Fuel type, HM mass, Power, time vector, and TMVA input variables names (looks the manual for format details)
	string fDBFType;		//!<  Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
	string fDBRType;		//!<  Reactor Type (e.g PWR, FBR-Na, ADS..)
	double 	fSpecificPower; 	//!< The specific power in W/gHM (HM: heavy Metal)
	double  fMaximalContent;	//!< The approx. maximum fissile content reachable by the model
	
	/*!
	 \name Others 
	 */
	//@{
	double fTotalFissileMassInStocks;//!< Total mass in fissile stocks
	double fTotalFertileMassInStocks;//!< Total mass in fertile stocks
	//@}
};

#endif









