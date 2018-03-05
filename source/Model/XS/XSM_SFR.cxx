#include <dirent.h>
#include <errno.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>


#include <TGraph.h>
#include <TString.h>
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "XSModel.hxx"
#include "XSM_SFR.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "StringLine.hxx"

//________________________________________________________________________
//
//        XSM_SFR
//
//________________________________________________________________________
XSM_SFR::XSM_SFR(string TMVA_Weight_Directory,
                 map<string, double> FixedParameters,
                 string InformationFile) {
    
    fLog = new CLASSLogger("XSM_SFR.log");
  
    DBGL
    fTMVAWeightFolder = TMVA_Weight_Directory;
    
    fInformationFile = fTMVAWeightFolder+InformationFile;
    
    DBGL
    GetMLPWeightFiles();
    DBGL
    
    INFO << "__A cross section interpolator using" << endl;
    INFO << "Multi Layer Perceptron has been define__" << endl;
    INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
    
    LoadKeyword();
    DBGL
    ReadNFO();
    DBGL
    BookTMVAReader(); 
    SetFixedVariablesValues(FixedParameters);
    DBGL
}
//________________________________________________________________________
XSM_SFR::XSM_SFR(CLASSLogger* Log,
                 string TMVA_Weight_Directory,
                 map<string,double> FixedParameters,
                 string InformationFile) {
    DBGL
    fLog = Log;
      
    fTMVAWeightFolder = TMVA_Weight_Directory;
    
    fInformationFile = fTMVAWeightFolder+InformationFile;
    
    DBGL
    GetMLPWeightFiles();
    DBGL
    
    INFO << "__A cross section interpolator using" << endl;
    INFO << "Multi Layer Perceptron has been define__" << endl;
    INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
    
    LoadKeyword();
    DBGL
    ReadNFO();
    DBGL
    BookTMVAReader(); 
    SetFixedVariablesValues(FixedParameters);
    DBGL
}



//________________________________________________________________________
XSM_SFR::~XSM_SFR()
{
    DBGL
    fTMVAVariableNames.clear();
    fDKeyword.clear();
    DBGL
}

//________________________________________________________________________
void XSM_SFR::SetFixedVariablesValues(map<string, double> FixedParameters) {
  
  fTMVAFixedVariableValues = vector<double>(fTMVAVariableNames.size(), -1);
  for (unsigned int i = 0; i < fTMVAVariableNames.size(); i++) {
    if (FixedParameters.find(fTMVAVariableNames[i]) != FixedParameters.end()) {
      fTMVAFixedVariableValues[i] = FixedParameters[fTMVAVariableNames[i]];
    }
  }
}

//________________________________________________________________________
void XSM_SFR::LoadKeyword()
{
    DBGL
    XSM_MLP::LoadKeyword();
    fDKeyword.insert( pair<string, XS_SFR_DMthPtr>( "k_timestep",    & XSM_MLP::ReadTimeSteps));
    fDKeyword.insert( pair<string, XS_SFR_DMthPtr>( "k_inputparameter",    & XSM_SFR::ReadNonZaiTMVAVariables)     );
    fDKeyword.insert( pair<string, XS_MLP_DMthPtr>( "k_zainame",    & XSM_MLP::ReadZAIName)     );
    DBGL
}
//________________________________________________________________________
void XSM_SFR::ReadNonZaiTMVAVariables(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_inputparameter" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_inputparameter\" not found !" << endl;
        exit(1);
    }
    
    string name = StringLine::NextWord(line, pos, ' ');
    
    fTMVAVariableNames.push_back( name );
    
    DBGL
}
//________________________________________________________________________
void XSM_SFR::ReadLine(string line)
{
    DBGL
    
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    
    map<string, XS_SFR_DMthPtr>::iterator it = fDKeyword.find(keyword);
    if(it != fDKeyword.end())
        (this->*(it->second))( line );
    
    DBGL
}
//________________________________________________________________________
void XSM_SFR::GetMLPWeightFiles()
{
    DBGL
    /**********Get All TMVA weight files*******************/
    //check if the folder containing weights exists
    DIR* rep = NULL;
    struct dirent* fichierLu = NULL;
    rep = opendir(fTMVAWeightFolder.c_str());
    
    if (rep ==  NULL)
    {
        ERROR << " Reading error for TMVA weight folder  " << fTMVAWeightFolder.c_str() << " : " << strerror(errno) << endl;
        exit(1);
    }
    
    /**Save file names of TMVA  weights*/
    fWeightFiles.resize(0);
    while ((fichierLu = readdir(rep)) != NULL)
    {
        string FileName = fichierLu->d_name ;
        if(FileName != "." && FileName != ".." )
        {
            if(FileName[FileName.size()-3] == 'x'  &&  FileName[FileName.size()-2] == 'm' && FileName[FileName.size()-1] == 'l' && FileName[0] != '.' )
                fWeightFiles.push_back(FileName);
            
        }
    }
    DBGL
}
//________________________________________________________________________
//________________________________________________________________________
//
//    Time  (MLP take time as parameter)
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

void XSM_SFR::BookTMVAReader() {
    DBGL

    string dir = fTMVAWeightFolder;
    if (dir[dir.size() - 1] != '/') {
        dir += "/";
    }

    for (int i = 0; i < int(fWeightFiles.size()); i++) {
        int Z = -2;
        int A = -2;
        int I = -2;
        int Reaction = -2;
        ReadWeightFile(fWeightFiles[i], Z, A, I, Reaction);
        if (Z >= GetZAIThreshold()) {
    DBGL
            fReader.push_back(new CLASSReader(fTMVAVariableNames));
    DBGL
            
            map<ZAI,string>::iterator it; 
            for ( it = fMapOfTMVAVariableNames.begin(); it != fMapOfTMVAVariableNames.end(); it++){
              fReader[i]->AddVariable(it->second.c_str());
            }
    DBGL
            // Time as to be the last one !!! 
            fReader[i]->AddVariable("Time");
    DBGL
            fReader[i]->BookMVA("MLP method", dir + fWeightFiles[i]);
    DBGL
        }
    }
    DBGL
}

vector<float> XSM_SFR::CreateTMVAInput(IsotopicVector IV, int TimeStep)
{
    DBGL
  vector<float> tmva_input_iv = XSM_MLP::CreateTMVAInput(IV, TimeStep);
    
  vector<float> tmva_input;  
  for(unsigned int j = 0 ; j < fTMVAVariableNames.size() ; j++) {
            tmva_input.push_back(fTMVAFixedVariableValues[j]);
    }    
  tmva_input.insert(tmva_input.end(), tmva_input_iv.begin(), tmva_input_iv.end());    
    return tmva_input;
}
//________________________________________________________________________
//________________________________________________________________________
