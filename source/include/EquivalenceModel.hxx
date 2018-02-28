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

    map < string, IsotopicVector> GetAllStreamList() {return fStreamList;}		// return all the lists

    virtual  map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListMassFractionMin, map < string , double> StreamListMassFractionMax, map < int , string> StreamListPriority, map < string , bool> StreamListIsBuffer);
    
	double SecondToBurnup(double Second) {return Second * fSpecificPower / (24 * 3.6e6);}
	double BurnupToSecond(double BurnUp) {return BurnUp / fSpecificPower * (24 * 3.6e6);}

	bool isIVInDomain(IsotopicVector IV);
	void StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks);
	void ConvertMassToLambdaVector(string MaterialDenomination, vector<double>& lambda, double MaterialMassNeeded, vector <IsotopicVector>  Stocks);

protected :

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









