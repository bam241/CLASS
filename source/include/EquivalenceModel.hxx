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
	EquivalenceModel(string TMVAXMLFilePath, string TMVANFOFilePath);			//!< Default constructor with path
	EquivalenceModel(CLASSLogger* log, string TMVAXMLFilePath, string TMVANFOFilePath);	//!< Logger constructor with path

	EquivalenceModel(string TMVANFOFilePath);			//!< Default constructor without Eq Model
	EquivalenceModel(CLASSLogger* log, string TMVANFOFilePath);	//!< Logger constructor Without Eq Model

	virtual ~EquivalenceModel();		//!< Destructor
	//@}

	string GetTMVAXMLFilePath() {return fTMVAXMLFilePath;} // Return the path to TMVA XML File path
	string GetTMVANFOFilePath() {return fTMVANFOFilePath;} // Return the path to TMVA NFO File path
	void SetTMVAXMLFilePath(string TMVAXMLFilePath) {fTMVAXMLFilePath = TMVAXMLFilePath;} // Set the path to TMVA XML File path
	void SetTMVANFOFilePath(string TMVANFOFilePath) {fTMVANFOFilePath = TMVANFOFilePath;} // Set the path to TMVA NFO File path

	double SecondToBurnup(double Second) {return Second * fSpecificPower / (24 * 3.6e6);}
	double BurnupToSecond(double BurnUp) {return BurnUp / fSpecificPower * (24 * 3.6e6);}

    bool isIVInDomain(IsotopicVector IV);
    void StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks);
    void ConvertMassToLambdaVector(string MaterialDenomination, vector<double>& lambda, double MaterialMassNeeded, vector <IsotopicVector>  Stocks);    

protected :

	bool fUseTMVAPredictor; //!< Bool that says if we need a TMVA predictor. If not, fuel fraction isimpoased by the FP.

	void SetLambdaToErrorCode(vector<double>& lambda);

	string fTMVAXMLFilePath;        //!<The weight needed by TMVA to construct and execute the multilayer perceptron
	string fTMVANFOFilePath;        //!<The weight needed by TMVA to construct and execute the multilayer perceptron

#ifndef __ROOTCLING__
	map<string, EQM_MthPtr> fKeyword;
#endif

	bool freaded;

	/*!
	 \name Others
	 */
	//@{
	map <string , double > fTotalMassInStocks;  		//!< Total mass in each vector of stock
	map <string , double > fLambdaMax;  			//!< Total lambda of available stocks

	//@}
};

#endif









