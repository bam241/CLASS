#ifndef _EQUIVALENCEMODEL_HXX
#define _EQUIVALENCEMODEL_HXX


/*!
 \file
 \brief Header file for EquivalenceModel class.
 
 
 @author BaM
 @author BLG
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
 @see EQM_BakerRoss_FBR_MOX
 @see EQM_LIN_PWR_MOX
 @see EQM_MLP_PWR_MOX
 @see EQM_POL_PWR_UO2
 @see EQM_QUAD_PWR_MOX
 
 @author BaM
 @author BLG
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
	 */
	virtual	 map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray);
	//}
		
	//@}
	
	/*!
	 \name Get/Set Method
	 */
	//@{

	IsotopicVector GetStreamList(string keyword) {return fStreamList[keyword];}	//!<return one of the list
	map < string, IsotopicVector> GetAllStreamList() {return fStreamList;}	 	//!<return allthe list

	int GetStreamListNumber(){return fStreamList.size();};
	
	map < string, double> GetBuildFuelFirstGuess(){return fFirstGuessContent;} 	//!<Get the initialization value for BuildFuel algorithm
	map < string, double> GetActualMassContent() const  { return  fActualMassContentInFuel ; }
	map < string, double> GetActualMolarContent() const { return fActualMolarContentInFuel ; }

	virtual map < string , double> GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp) = 0; //!< Return the molar fraction of fissile element in the fuel according to the burnup, and a given fuel composition (this is the heart of the equivalence model)

	double GetRelativMassPrecision() const	{ return fRelativMassPrecision; }	//!< Mass precision
	int GetMaxInterration()		 const	{ return fMaxInterration; }		//!< Max iterration in build fueld algorythm
	
	void SetBuildFuelFirstGuess(map < string, double> FirstGuess){fFirstGuessContent = FirstGuess;} //!<set the initialization value for BuildFuel algorithm
	void SetRelativMassPrecision( double val)	{ fRelativMassPrecision = val; }	//!< Mass precision
	void SetMaxInterration(int val)			{ fMaxInterration = val; }	//!< Max iterration in build fueld algorythm

	
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
	
	//@}
	
	
	bool isIVInDomain(IsotopicVector IV);

	
	
	
	protected :
	
	map < string, IsotopicVector> fStreamList;		//!< contains all lists of zai needed to build a fuel (example : 2 -> fissileList+fertileList)
								//!< each list is identified by a keyword (example : -> "Fissil" & "Fertil")
								
	
	map < string, double> fFirstGuessContent;		//!< fissile content for BuildFuel initialization (in weight proportion)
	
	map < string, double> fActualMolarContentInFuel; 	//!< Molar Content in fuel of each list at this step of the calculation
	map < string, double> fActualMassContentInFuel; 	//!< Molar Content in fuel of each list at this step of the calculation

	
	double LambdaCalculation(string MaterialDenomination, double MaterialMassNeeded, vector <IsotopicVector>  Stocks); //!< Calculate the proportion of each stocks in StockArray to take in oder to get a mass of MassNeeded (can be Fer(fertile) or Fis(Fissile))
	void SetLambda(vector<double>& lambda , double Lambda_tot);	//!< Set individual lambda according to the LAMBDA_TOT (lambda of all stocks)	
	void SetLambdaToErrorCode(vector<double>& lambda);
	void StocksTotalMassCalculation(map < string , vector <IsotopicVector> >  Stocks);
	double fRelativMassPrecision;				//!< Mass precision
	int fMaxInterration;					//!< Max iterration in build fueld algorythm

	
	//@}

#ifndef __CINT__
	map<string, EQM_MthPtr> fKeyword;
#endif
	
	bool freaded;
	map< ZAI, pair<double,double> > fZAILimits; //!< Fresh fuel range : map<ZAI<min edge ,max edge >>

	string fInformationFile;		//!<  file containing Reactor Type, Fuel type, HM mass, Power, time vector, and TMVA input variables names (looks the manual for format details)
	string fDBFType;		//!<  Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
	string fDBRType;		//!<  Reactor Type (e.g PWR, FBR-Na, ADS..)
	double 	fSpecificPower; 	//!< The specific power in W/gHM (HM: heavy Metal)
	double  fMaximalContent;	//!< The approx. maximum fissile content reachable by the model
	
	/*!
	 \name Others 
	 */
	//@{
	map <string , double > fTotalMassInStocks;  //!< Total mass in each vector of stock
	//@}



};

#endif









