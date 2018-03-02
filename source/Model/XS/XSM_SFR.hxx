#ifndef _XSM_SFR_HXX
#define _XSM_SFR_HXX

/*!
 \file
 \brief Header file for XSM_SFR class.
 
 
 @authors Marc
 @version 1.0
 */
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "TTree.h"

#include "XSM_MLP.hxx"

typedef long long int cSecond;
using namespace std;


class XSM_SFR;
#ifndef __CINT__
typedef void (XSM_SFR::*XS_SFR_DMthPtr)( const string & ) ;
#endif
//-----------------------------------------------------------------------------//
//! Defines a XSModel getting mean cross sections from neural network execution

/*!
 Define a XSM_SFR.
 This is the class to predict cross sections with a
 set of Multi Layer Perceptrons (MLP)
 
 @authors Marc
 @version 1.0
 */
//________________________________________________________________________


class XSM_SFR : public XSM_MLP
{
    public :
    
    /*!
     \name Constructor/Desctructor
     */
    //@{
    
    //{
    /// Normal Constructor
    /*!
     \param TMVA_Weight_Directory : The directory where all the TMVA weight are located
     \param InformationFile : Name of the information file located in TMVA_Weight_Directory (default : Data_Base_Info.nfo)
     \param IsTimeStep : if true , one TMVA weihgt per step time is requiered otherwise it assumes time is part of the MLP inputs
     
     */
    XSM_SFR(string TMVA_Weight_Directory,map<string,double> FixedParameters,string InformationFile = "/Data_Base_Info.nfo");
    //}
    
    //{
    /// CLASSLogger Constructor
    /*!
     \param log : The CLASSLogger
     \param TMVA_Weight_Directory : The directory where all the TMVA weight are located
     \param InformationFile : Name of the information file located in TMVA_Weight_Directory (default : Data_Base_Info.nfo)
     \param IsTimeStep : if true , one TMVA weihgt per step time is requiered otherwise it assumes time is part of the MLP inputs
     
     */
    XSM_SFR(CLASSLogger* Log,string TMVA_Weight_Directory,map<string,double> FixedParameters,string InformationFile = "/Data_Base_Info.nfo");
    //}
    
    ~XSM_SFR();
    //@}
    
    /*!
     \name Reading NFO related Method
     */
    //@{
    
    //{
    /// LoadKeyword() : make the correspondance between keyword and reading method
    void LoadKeyword();
    //}
    
    //{
    /// ReadTimeSteps : read the time step of the model
    /*!
     \param line : line suppossed to contain the time step information starts with "k_timestep" keyword
     */
    void ReadTimeSteps(const string &line);
    //}
    
    //{
    /// ReadZAIName : read the zai name in the TMWA MLP model
    /*!
     \param line : line suppossed to contain the ZAI name  starts with "k_zainame" keyword
     */
    void ReadNonZaiTMVAVariables(const string &line);
    //}
    //{
    /// ReadLine : read a line
    /*!
     \param line : line to read
     */
    void ReadLine(string line);
    //}
    
    //@}
    void FixTMVAVariable(string VariableName,double VariableValue);

    void SetFixedVariablesValues(map<string,double> FixedParameters);
    void BookTMVAReader();
    
    private :
    
    void GetMLPWeightFiles();                //!< Find all .xml file in TMVA_Weight_Directory
    vector<float> CreateTMVAInput(IsotopicVector isotopicvector,int TimeStep = 0);    //!<Create input tmva to be read by ExecuteTMVA
    
    vector<double>     fMLP_Time;    //!<  Time vector of the data base
    vector<string>     fWeightFiles;    //!<  All the weight file contains in fTMVAWeightFolder
    string fTMVAWeightFolder;    //!<  folder containing all the weight file
    bool fIsStepTime;        //!<  true if one TMVA weihgt per step time is requiered otherwise it assumes time is part of the MLP inputs
    vector<string> fTMVAVariableNames;//!<  List of TMVA input variable names (read from fMLPInformationFile ) , name depends on the training step
    vector<bool> fTMVAFixedVariable;//!< List of value for TMVA that have to be used for all 
    vector<double> fTMVAFixedVariableValues;//!< List of value for TMVA that have to be used for all 
    
#ifndef __CINT__
    map<string, XS_SFR_DMthPtr> fDKeyword;
#endif
};

#endif

