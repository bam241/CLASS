#ifndef _EQUIVALENCE_MODEL_
#define _EQUIVALENCE_MODEL_


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
    EquivalenceModel();            //!< Default constructor with path
    EquivalenceModel(CLASSLogger* log);    //!< Logger constructor with path

    virtual ~EquivalenceModel();        //!< Destructor
    //@}

    map < string, IsotopicVector> GetAllStreamList() {return fStreamList;}        // return all the lists

    virtual  map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListMassFractionMin, map < string , double> StreamListMassFractionMax, map < int , string> StreamListPriority, map < string , bool> StreamListIsBuffer); //!< Main function to build the fuel of mass HMMAss that will reach BurnUp in the reactor -> Function that have to be defined in source/Model/Equivalence
    
    double SecondToBurnup(double Second) {return Second * fSpecificPower / (24 * 3.6e6);}
    double BurnupToSecond(double BurnUp) {return BurnUp / fSpecificPower * (24 * 3.6e6);}



    void LoadKeyword();
    void ReadNFO();
      virtual void ReadLine(string line);
    void ReadZAIlimits(const string &line);
    void ReadType(const string &line);
    
    void ReadList(const string &line);

    void ReadEqMaxFraction(const string &line); 
    void ReadEqMinFraction(const string &line);
    //{
    /// ReadZAIName : read the zai name in the TMWA MLP model
    /*!
     \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
     */
    void ReadZAIName(const string &line);
    //}
    
    //{
    /// ReadPredictorType: read the type of predictor used (ex : MLP)
    /*!
     \param line : line suppossed to contain the Buffer information starts with "k_predictortype" keyword
     */
    void ReadPredictorType(const string &line);
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

    //{
    /// ReadTargetParameter : type of target parameter optimized in build fuel (ex. BUmax)
    /*!
     \param line : line suppossed to contain the Target Parameter information starts with "k_targetparameter" keyword
     */
    void ReadTargetParameter(const string &line);
    //}
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

    //{
    /// ReadBuffer : read the Buffer material name in the fuel
    /*!
     \param line : line suppossed to contain the Buffer information starts with "k_buffer" keyword
     */
    void ReadBuffer(const string &line);
    //}




    bool isIVInDomain(IsotopicVector IV); //<!< Check if the current IV is in the domain allowed for the fuel fabrication
    void StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks); //!< Calculate the Total Mass in the dStock 
    void ConvertMassToLambdaVector(string MaterialDenomination, vector<double>& lambda, double MaterialMassNeeded, vector <IsotopicVector>  Stocks); //!< Calculate the fraction of IV that should be take to build the fuel (Lambda = 1 <=> the whole IV is used for the fuel  

protected :

    map < string, IsotopicVector> fStreamList;                  //!< contains all lists of zai needed to build a fuel (example : 2 -> fissileList+fertileList)
                                                                //!< each list is identified by a keyword (example : -> "Fissile" & "Fertile")
    double     fSpecificPower;                                      //!< The specific power in W/gHM (HM: heavy Metal)
    void SetLambdaToErrorCode(vector<double>& lambda);


#ifndef __ROOTCLING__
    map<string, EQM_MthPtr> fKeyword;    //!< Parameters of the equivalence model 
#endif

    bool freaded;

    map< ZAI, pair<double, double> > fZAILimits;     //!< Fresh fuel range : map<ZAI<min edge ,max edge >>

    /*!
     \name Others
     */
    //@{
    map <string , double > fTotalMassInStocks;      //!< Total mass in each vector of stock
    map <string , double > fLambdaMax;              //!< Total lambda of available stocks

    //@}


    string fInformationFile;                    //!<  file containing Reactor Type, Fuel type, HM mass, Power, time vector, and TMVA input variables names (looks the manual for format details)
    string fDBFType;                    //!<  Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
    string fDBRType;                    //!<  Reactor Type (e.g PWR, FBR-Na, ADS..)

    
    map<ZAI, string> fMapOfTMVAVariableNames;               //!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step

    double  fMaximalBU;                                 //!< The Maximum burn-up of the model in GWd/t
    map < string , double> fStreamListEqMMassFractionMax;           //!< Map that contains lists of stream according to the EqModel with mass maximum fraction
    map < string , double> fStreamListEqMMassFractionMin;           //!< Map that contains lists of stream according to the EqModel with mass minimum fraction

    string fPredictorType ;                                 //!< Type of predictor used in Equivalence Model (ex: MLP)
    string fOutput ;                                //!< Type of output calculated by the predictor
    string fBuffer ;                                    //!< Name of material used as buffer in fuel
};

#endif









