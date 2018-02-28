#ifndef _EQONEPARAMETER_HXX_
#define _EQONEPARAMETER_HXX_

/*!
 \file
 \brief Header file for EQ_OneParameter class.
 
 
 @author BLG
 @author BaM
 @author FaC
 @version 3.0
 */

#include "EquivalenceModel.hxx"
#include "IsotopicVector.hxx"
#include "CLASSReader.hxx"
#include <math.h>
#include "TTree.h"
#include <map>


using namespace std;

class EQ_OneParameter;
typedef void (EQ_OneParameter::*EQOP_MthPtr)( const string & );

//-----------------------------------------------------------------------------//

//! Determines how to build a fresh fuel 
/*!
 Define an EQ_OneParameter.
 The aim of these class is to gather all the commum properties of all
 Equivalence Model.
 
 \warning
 Never instantiate EQ_OneParameter in your CLASS input but it's derivated class
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


class EQ_OneParameter : public EquivalenceModel
{
    public :
    /*!
     \name Constructor/Desctructor
     */
    //@{
    EQ_OneParameter(string TMVAXMLFilePath, string TMVANFOFilePath);           //!< Default constructor with path
    EQ_OneParameter(CLASSLogger* log, string TMVAXMLFilePath, string TMVANFOFilePath); //!< Logger constructor with path
    
    EQ_OneParameter(string TMVANFOFilePath);           //!< Default constructor without Eq Model
    EQ_OneParameter(CLASSLogger* log, string TMVANFOFilePath); //!< Logger constructor Without Eq Model

    virtual ~EQ_OneParameter();        //!< Destructor
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
    virtual  map <string , vector<double> > BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListMassFractionMin, map < string , double> StreamListMassFractionMax, map < int , string> StreamListPriority, map < string , bool> StreamListIsBuffer);
    //}


    //@}
    void BookTMVAReader(); //Book TMVA method  
    
    vector<float> CreateTMVAInput(IsotopicVector TheFreshfuel, double ThisTime)    ; //!<Create input tmva tree to be read by ExecuteTMVA
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

    int GetStreamListNumber(){return fStreamList.size();};
    int GetMaxIterration()       const  { return fMaxIterration; }      //!< Max iterration in build fueld algorythm    
    double GetTargetParameterStDev(){return fTargetParameterStDev;}//!< Get the precision on fTargetParameterStDev
    double GetStreamListEqMMassFractionMax(string keyword){return fStreamListEqMMassFractionMax[keyword] ;}
    double GetStreamListEqMMassFractionMin(string keyword){return fStreamListEqMMassFractionMin[keyword] ;}
    double GetPCMPrecision(){return fPCMprecision/1e5;}//!< Get the precision on @f$\langle k \rangle@f$ prediction []. Neural network predictor constructors

    void SetModelParameter(string sMP, double dMP)  { fModelParameter[sMP] = dMP; }   //!< Set Model Parameters precised in NFO file
    map<string, double> GetModelParameter()  { return fModelParameter; }   //!< Get Model Parameters precised in NFO file
	  
    void SetNonZaiTMVAVariable(string snZP, double dnZP);    //!< Set NonZaiTMVAVariables
	vector< pair<double, string> > GetNonZaiTMVAVariables()  { return fListOfNonZaiTMVAVariables; }   //!< Get NonZaiTMVAVariables

    void SetMaxIterration(int val)  { fMaxIterration = val; }   //!< Max iteration in build fuel algorithm
    void SetTargetParameterStDev(double TPSD){fTargetParameterStDev = TPSD;} //!< Set the precision on Target Parameter
    void SetStreamListEqMMassFractionMax(string keyword, double value){fStreamListEqMMassFractionMax[keyword] = value;}
    void SetStreamListEqMMassFractionMin(string keyword, double value){fStreamListEqMMassFractionMin[keyword] = value;}

    void SetPCMPrecision(double pcm){fPCMprecision = pcm;}        //!< Set the precision on @f$\langle k \rangle@f$ prediction [pcm]. Neural network predictor constructors

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

    //{
    /// ReadModelParameter : read the name of equivalence model parameter
    /*!
     \param line : line suppossed to contain the Buffer information starts with "k_modelparameter" keyword
     */
    void ReadModelParameter(const string &line);
    //} 

    //{
    /// ReadPredictorType: read the type of predictor used (ex : MLP)
    /*!
     \param line : line suppossed to contain the Buffer information starts with "k_predictortype" keyword
     */
    void ReadPredictorType(const string &line);
    //}

    //{
    /// ReadTargetParameterStDev: read the target parameter standard deviation
    /*!
     \param line : line suppossed to contain the Buffer information starts with "k_targetparameterstdev" keyword
     */
    void ReadTargetParameterStDev(const string &line);
    //}

    void PrintInfo(); //Print the information red in the INFO stream    
    
    //{
    /// ReadFissil : read the zai name in the TMWA MLP model starts with "k_fissil" keyword
    /*!
     \param line : line suppossed to contain the fissil list
     */
    void ReadList(const string &line);

    void ReadEqMaxFraction(const string &line); 
    void ReadEqMinFraction(const string &line);
    
    IsotopicVector BuildFuelToTest(map < string, vector<double> >& lambda, map < string , vector <IsotopicVector> > const& StreamArray, double HMMass, map <string, bool> StreamListIsBuffer); //Build a fuel with the buffer according to fissile lambda
    void CheckTargetParameterConsistency(map < int , string >  StreamListPriority, map < string , double >  TargetParameterMin, map < string , double > TargetParameterMax);

    protected :

    bool fUseTMVAPredictor; //!< Bool that says if we need a TMVA predictor. If not, fuel fraction isimpoased by the FP.
    
    map < string , double> fStreamListEqMMassFractionMax;           //!< Map that contains lists of stream according to the EqModel with mass maximum fraction
    map < string , double> fStreamListEqMMassFractionMin;           //!< Map that contains lists of stream according to the EqModel with mass minimum fraction

    string fPredictorType ;                                 //!< Type of predictor used in Equivalence Model (ex: MLP)
    string fOutput ;                                //!< Type of output calculated by the predictor
    string fBuffer ;                                    //!< Name of material used as buffer in fuel

    map<string, double> fModelParameter ;                   ///< Map of equivalence model parameter
 	  vector< pair<double, string> > fListOfNonZaiTMVAVariables ; 					///!<  List of TMVA input variable names that are not ZAIs

    map<ZAI, string> fMapOfTMVAVariableNames;               //!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step

    double  fTargetParameterStDev;                          //!< Precision on target parameter calculation

    double  fMaximalBU;                                 //!< The Maximum burn-up of the model in GWd/t
    string fTargetParameter;                            //!< Type of target parameter optimized in build fuel (ex. BUmax)               
    int fMaxIterration;                             //!< Max iterrations in build fueld algorithm    

    string fTMVAXMLFilePath;        //!<The weight needed by TMVA to construct and execute the multilayer perceptron
    string fTMVANFOFilePath;        //!<The weight needed by TMVA to construct and execute the multilayer perceptron


    CLASSReader* fReader;
    map<string, EQOP_MthPtr> fKeyword;
    
    map< ZAI, pair<double,double> > fZAILimits;     //!< Fresh fuel range : map<ZAI<min edge ,max edge >>

    string fInformationFile;                    //!<  file containing Reactor Type, Fuel type, HM mass, Power, time vector, and TMVA input variables names (looks the manual for format details)
    string fDBFType;                    //!<  Fuel Type    (e.g MOX, UOX, ThU, ThPu ...)
    string fDBRType;                    //!<  Reactor Type (e.g PWR, FBR-Na, ADS..)

    private :

    double fPCMprecision;          //!< precision on @f$\langle k \rangle@f$ prediction [pcm]

};

#endif









