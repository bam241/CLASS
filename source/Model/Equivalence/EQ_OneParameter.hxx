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

#include <math.h>
#include <map>
#include "CLASSReader.hxx"
#include "EquivalenceModel.hxx"
#include "IsotopicVector.hxx"
#include "TTree.h"

using namespace std;

class EQ_OneParameter;
typedef void (EQ_OneParameter::*EQ_OP_MthPtr)(const string&);

//-----------------------------------------------------------------------------//

//! Determines how to build a fresh fuel
/*!
 Define an EQ_OneParameter.
 The aim of this class is to build a fuel based on the optimisation of one
 parameter
 This one may be
     - keff at begining of cycle
    - The maximum achievable Burn Up
 For those two possibilities one may define some external parameters like th
 kthreshold of the number of batch for instance

 In case one wants to build a fuel with a constant fraction of fissile, it has
 to be set in the fabrication plant


 @author BLG
 @author BaM
 @author FaC

 @version 3.0
 */
//________________________________________________________________________

class EQ_OneParameter : public EquivalenceModel {
 public:
  /*!
   \name Constructor/Desctructor
   */
  //@{
  EQ_OneParameter(
      string TMVAXMLFilePath,
      string TMVANFOFilePath);  //!< Default constructor with path of TMVA files
  EQ_OneParameter(CLASSLogger* log, string TMVAXMLFilePath,
                  string TMVANFOFilePath);  //!< Logger constructor with path

  EQ_OneParameter(
      string TMVANFOFilePath);  //!< Default constructor without Eq Model
  EQ_OneParameter(
      CLASSLogger* log,
      string TMVANFOFilePath);  //!< Logger constructor Without Eq Model

  virtual ~EQ_OneParameter();  //!< Destructor
  //@}

  /*!
   \name Fuel Construction Method
   */
  //@{

  //{
  /// BuildFuel function.
  /*!
   Build the fuel following the equivalance model with the proper requierment in
   term of mass, burnup....
   \param double burnup reached by the fuel at the end of irradiation
   \param double HMMass, Heavy metal mass needed
   \param map < string , vector <IsotopicVector> > StreamArray, the string is
   the stream code (fissile fertile ,...) the IsotopicVector the fraction of
   each IV to take in the (fissile, fertile,..) stock .
   \param map < string , double> StreamListMassFractionMin defines the lower
   limits allowed for each Stream (in terms of percentage)
   \param map < string , double> StreamListMassFractionMax defines the upper
   limits allowed for each Stream
   \param map < int , string> StreamListPriority defines the priority of
   different streams (1 is the higher priority, then 2, then 3,... Two different
   streams can not have the same priority
   \param map < string , bool> StreamListIsBuffer defines if each stream is the
   fuel buffer (only one stream can be the fuel buffer)
   */
  virtual map<string, vector<double> > BuildFuel(
      double BurnUp, double HMMass,
      map<string, vector<IsotopicVector> > StreamArray,
      map<string, double> StreamListMassFractionMin,
      map<string, double> StreamListMassFractionMax,
      map<int, string> StreamListPriority,
      map<string, bool> StreamListIsBuffer);
  //}

  //@}
  void BookTMVAReader();  // Book TMVA method

  vector<float> CreateTMVAInput(
      IsotopicVector TheFreshfuel,
      double ThisTime);  //!<Create input tmva tree to be read by ExecuteTMVA
  double CalculateTargetParameter(IsotopicVector TheFuel,
                                  string TargetParameterName);  //!<Get a fuel
                                                                //! parameter
  //! associated to
  //! the fuel --->
  double CalculateBurnUpMax(IsotopicVector TheFuel);  //!<Calculate
                                                                  //! the BU max
  //! associated
  //! to a fuel
  //! composition
  //! based on
  //! MLP
  //! prediction
  //!(suitable
  //! for PWR)
  double CalculateKeffAtBOC(IsotopicVector TheFuel);  //!<Calculate the keff at
                                                      //! BOC associated to a
  //! fuel composition based
  //! on MLP prediction
  //!(suitable for SFR)

  /*!
   \name Get/Set Method
   */
  //@{

  int GetStreamListNumber() { return fStreamList.size(); };
  int GetMaxIterration() const {
    return fMaxIterration;
  }  //!< Max iterration in build fueld algorythm
  double GetTargetParameterStDev() {
    return fTargetParameterStDev;
  }  //!< Get the precision on fTargetParameterStDev
  double GetStreamListEqMMassFractionMax(string keyword) {
    return fStreamListEqMMassFractionMax[keyword];
  }
  double GetStreamListEqMMassFractionMin(string keyword) {
    return fStreamListEqMMassFractionMin[keyword];
  }
  double GetPCMPrecision() {
    return fPCMprecision / 1e5;
  }  //!< Get the precision on @f$\langle k \rangle@f$ prediction []. Neural
     //! network predictor constructors

  void SetModelParameter(string sMP, double dMP) {
    fModelParameter[sMP] = dMP;
  }  //!< Set Model Parameters precised in NFO file can be keff or BU Max
  map<string, double> GetModelParameter() {
    return fModelParameter;
  }  //!< Get Model Parameters precised in NFO file

  void SetNonZaiTMVAVariable(
      string snZP,
      double dnZP);  //!< Set NonZaiTMVAVariables for the input of the TMVA
  vector<pair<double, string> > GetNonZaiTMVAVariables() {
    return fListOfNonZaiTMVAVariables;
  }  //!< Get NonZaiTMVAVariables

  void SetMaxIterration(int val) {
    fMaxIterration = val;
  }  //!< Max iteration in build fuel algorithm
  void SetTargetParameterStDev(double TPSD) {
    fTargetParameterStDev = TPSD;
  }  //!< Set the precision on Target Parameter
  void SetStreamListEqMMassFractionMax(string keyword, double value) {
    fStreamListEqMMassFractionMax[keyword] = value;
  }
  void SetStreamListEqMMassFractionMin(string keyword, double value) {
    fStreamListEqMMassFractionMin[keyword] = value;
  }

  void SetPCMPrecision(double pcm) {
    fPCMprecision = pcm;
  }  //!< Set the precision on @f$\langle k \rangle@f$ prediction [pcm]. Neural
     //! network predictor constructors

  /*!
   \name Time <-> Burnup conversion
   */
  //@{

  double SecondToBurnup(double Second) {
    return Second * fSpecificPower / (24 * 3.6e6);
  }
  double BurnupToSecond(double BurnUp) {
    return BurnUp / fSpecificPower * (24 * 3.6e6);
  }

  //@}

  //@}

  /*!
   \name Reading NFO related Method
   */
  //@{
  void ReadLine(string line);
  void LoadKeyword();

  //{
  /// ReadModelParameter : read the name of equivalence model parameter
  /*!
   \param line : line suppossed to contain the Buffer information starts with
   "k_modelparameter" keyword
   */
  void ReadModelParameter(const string& line);
  //}

  //{
  /// ReadTargetParameter : type of target parameter optimized in build fuel (ex. BUmax)
  /*!
   \param line : line suppossed to contain the Target Parameter information
   starts with "k_targetparameter" keyword
   */
  void ReadTargetParameter(const string& line);
  //}

  //{
  /// ReadTargetParameterStDev: read the target parameter standard deviation
  /*!
   \param line : line suppossed to contain the Buffer information starts with
   "k_targetparameterstdev" keyword
   */
  void ReadTargetParameterStDev(const string& line);
  //}

  //{
  /// ReadNonZaiTMVAVariables : read the NonZai variables for the predictor (ex : Nbatch, Specific power)
  /*!
   \param line : line suppossed to contain the NonZai variables for TMVA starts
   with "k_nonZAIforTMVA" keyword
   */
  void ReadNonZaiTMVAVariables(const string& line);
  //}
  void PrintInfo();  // Print the information red in the INFO stream

  void CheckParameterValidity();

  IsotopicVector BuildFuelToTest(
      map<string, vector<double> >& lambda,
      map<string, vector<IsotopicVector> > const& StreamArray, double HMMass,
      map<string, bool> StreamListIsBuffer);  // Build a fuel with the buffer
                                              // according to fissile lambda
  void CheckTargetParameterConsistency(map<int, string> StreamListPriority,
                                       map<string, double> TargetParameterMin,
                                       map<string, double> TargetParameterMax);

 protected:
  bool fUseTMVAPredictor;  //!< Bool that says if we need a TMVA predictor. If
                           //! not, fuel fraction isimpoased by the FP.

  map<string, double> fModelParameter;  ///< Map of equivalence model parameter
  vector<pair<double, string> >
      fListOfNonZaiTMVAVariables;  ///!<  List of TMVA input variable names that are not ZAIs

  double fTargetParameterStDev;  //!< Precision on target parameter calculation

  string fTargetParameter;  //!< Type of target parameter optimized in build
                            //! fuel (ex. BUmax)
  int fMaxIterration;       //!< Max iterrations in build fueld algorithm

  CLASSReader* fReader;
  map<string, EQ_OP_MthPtr> fKeyword;  //!< The model parameters

 private:
  double
      fPCMprecision;  //!< precision on @f$\langle k \rangle@f$ prediction [pcm]
};

#endif
