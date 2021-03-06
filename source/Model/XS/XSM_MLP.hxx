#ifndef _XSM_MLP_HXX
#define _XSM_MLP_HXX

/*!
 \file
 \brief Header file for XSM_MLP class.


 @authors BaM
 @authors BLG
 @version 1.0
 */
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "CLASSReader.hxx"
#include "TTree.h"
#include "XSModel.hxx"

typedef long long int cSecond;
using namespace std;

class XSM_MLP;
#ifndef __ROOTCLING__
typedef void (XSM_MLP::*XS_MLP_DMthPtr)(const string &);
#endif
//-----------------------------------------------------------------------------//
//! Defines a XSModel getting mean cross sections from neural network execution

/*!
 Define a XSM_MLP.
 This is the class to predict cross sections with a
 set of Multi Layer Perceptrons (MLP)

 @authors BLG
 @version 1.0
 */
//________________________________________________________________________

class XSM_MLP : public XSModel {
 public:
  /*!
   \name Constructor/Desctructor
   */
  //@{
  XSM_MLP(){};

  //{
  /// Normal Constructor
  /*!
   \param TMVA_Weight_Directory : The directory where all the TMVA weight are
   located
   \param InformationFile : Name of the information file located in
   TMVA_Weight_Directory (default : Data_Base_Info.nfo)

   */
  XSM_MLP(string TMVA_Weight_Directory,
          string InformationFile = "/Data_Base_Info.nfo");
  //}

  //{
  /// CLASSLogger Constructor
  /*!
   \param log : The CLASSLogger
   \param TMVA_Weight_Directory : The directory where all the TMVA weight are
   located
   \param InformationFile : Name of the information file located in
   TMVA_Weight_Directory (default : Data_Base_Info.nfo)

   */
  XSM_MLP(CLASSLogger *Log, string TMVA_Weight_Directory,
          string InformationFile = "/Data_Base_Info.nfo");
  //}

  ~XSM_MLP();
  //@}

  /*!
   \name Reading NFO related Method
   */
  //@{

  //{
  /// LoadKeyword() : make the correspondance between keyword and reading method
  virtual void LoadKeyword();
  //}

  //{
  /// ReadTimeSteps : read the time step of the model
  /*!
   \param line : line suppossed to contain the time step information starts with
   "k_timestep" keyword
   */
  void ReadTimeSteps(const string &line);
  //}

  //{
  /// ReadZAIName : read the zai name in the TMWA MLP model
  /*!
   \param line : line suppossed to contain the ZAI name  starts with "k_zainame"
   keyword
   */
  void ReadZAIName(const string &line);
  //}
  //{
  /// ReadLine : read a line
  /*!
   \param line : line to read
   */
  void ReadLine(string line);
  //}

  //@}

  virtual void BookTMVAReader();  // Book TMVA method

  EvolutionData GetCrossSections(
      IsotopicVector IV,
      double t = 0);  //!< Return calculated cross section by the MLP regression

 protected:
  vector<CLASSReader *> fReader;

  void GetMLPWeightFiles();  //!< Find all .xml file in TMVA_Weight_Directory

  void ReadWeightFile(string Filename, int &Z, int &A, int &I,
                      int &Reaction);  //!<  Select the reaction according to
                                       //! the weight file name

  double ExecuteTMVA(string WeightFile, TTree *InputTree);  //!<Execute the MLP
                                                            //! according to the
  //! input tree
  //! created
  virtual vector<float> CreateTMVAInput(
      IsotopicVector isotopicvector,
      int t);  //!<Create input tmva tree to be read by ExecuteTMVA

  vector<double> fMLP_Time;  //!<  Time vector of the data base
  vector<string>
      fWeightFiles;  //!<  All the weight file contains in fTMVAWeightFolder

  string fTMVAWeightFolder;  //!<  folder containing all the weight file

  bool fIsStepTime;  //!<  true if one TMVA weihgt per step time is requiered
                     //! otherwise it assumes time is part of the MLP inputs

  map<ZAI, string> fMapOfTMVAVariableNames;  //!<  List of TMVA input variable
                                             //! names (read from
//! fMLPInformationFile ) , name
//! depends on the training step

#ifndef __ROOTCLING__
  map<string, XS_MLP_DMthPtr> fDKeyword;
#endif
};

#endif
