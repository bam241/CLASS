
#include "XSModel.hxx"
#include "XSM_MLP.hxx"
#include "CLASSLogger.hxx"
#include "CLASSMethod.hxx"
#include "CLASSReader.hxx"
#include "StringLine.hxx"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <TGraph.h>
#include <TString.h>
#include "TFile.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include <dirent.h>
#include <errno.h>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

//________________________________________________________________________
//
//        XSM_MLP
//
//
//
//
//________________________________________________________________________
XSM_MLP::XSM_MLP(string TMVA_Weight_Directory,string InformationFile):XSModel(new CLASSLogger("XSM_MLP.log"))
{
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
    DBGL
    
}
//________________________________________________________________________
XSM_MLP::XSM_MLP(CLASSLogger* Log,string TMVA_Weight_Directory,string InformationFile):XSModel(Log)
{
    DBGL
    fTMVAWeightFolder = TMVA_Weight_Directory;

    fInformationFile = fTMVAWeightFolder+InformationFile;
    
    GetMLPWeightFiles();
    DBGL
    
    INFO << "__A cross section interpolator using" << endl;
    INFO << "Multi Layer Perceptron has been define__" << endl;
    INFO << " \t His TMVA folder is : \" " << fTMVAWeightFolder << "\"" << endl;
    
    DBGL
    LoadKeyword();
    DBGL
    ReadNFO();
    DBGL
    BookTMVAReader(); 
    DBGL
}
//________________________________________________________________________
XSM_MLP::~XSM_MLP()
{
    DBGL
    fMapOfTMVAVariableNames.clear();
    fDKeyword.clear();
    //delete fReader;
    DBGL
}


//________________________________________________________________________
void XSM_MLP::LoadKeyword()
{
    DBGL
    XSModel::LoadKeyword();
    fDKeyword.insert( pair<string, XS_MLP_DMthPtr>( "k_timestep",    & XSM_MLP::ReadTimeSteps));
    fDKeyword.insert( pair<string, XS_MLP_DMthPtr>( "k_zainame",    & XSM_MLP::ReadZAIName)     );
    DBGL
}


//________________________________________________________________________
void XSM_MLP::ReadZAIName(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_zainame" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_zainame\" not found !" << endl;
        exit(1);
    }
    
    int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int I = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    
    string name = StringLine::NextWord(line, pos, ' ');
    
    fMapOfTMVAVariableNames.insert( pair<ZAI,string>( ZAI(Z, A, I), name ) );
    
    DBGL
}

//________________________________________________________________________
void XSM_MLP::ReadTimeSteps(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_timestep" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_timestep\" not found !" << endl;
        exit(1);
    }
    
    while( pos < (int)line.size() )
        fMLP_Time.push_back( atof( (StringLine::NextWord(line,pos,' ')).c_str() ));
    DBGL
}

//________________________________________________________________________
void XSM_MLP::ReadLine(string line)
{
    DBGL
    
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    
    map<string, XS_MLP_DMthPtr>::iterator it = fDKeyword.find(keyword);
    
    if(it != fDKeyword.end())
        (this->*(it->second))( line );
    
    DBGL
}


 
//________________________________________________________________________
void XSM_MLP::GetMLPWeightFiles()
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
void XSM_MLP::ReadWeightFile(string Filename, int &Z, int &A, int &I, int &Reaction)
{
    Z = -1;
    A = -1;
    I = -1;
    Reaction = -1;
    
    size_t found = Filename.find("XS");
    
    string NameJOB;
    NameJOB = Filename.substr(found);
    int pos = 0;
    
    StringLine::NextWord(NameJOB, pos, '_');
    
    Z = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
    
    A = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
    
    I = atof( (StringLine::NextWord(NameJOB,pos,'_') ).c_str() );
    
    string  sReaction = (StringLine::NextWord(NameJOB,pos,'_') ).c_str() ;
    size_t foundext = sReaction.find(".weights.xml");
    sReaction = sReaction.substr(0,foundext);
    
    if(sReaction == "fis")
        Reaction = 0;
    if(sReaction == "cap")
        Reaction = 1;
    if(sReaction == "n2n")
        Reaction = 2;
    
    if(Z<= 0 || A<= 0 || I<0 || Reaction == -1)
    {
        ERROR << " wrong TMVA weight format " << endl;
        exit(0);
    }
    
}

void XSM_MLP::BookTMVAReader() {

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
            fReader.push_back(new CLASSReader(fMapOfTMVAVariableNames));
            fReader[i]->AddVariable("Time");
            fReader[i]->BookMVA("MLP method", dir + fWeightFiles[i]);
        }
    }
    DBGL
}

//________________________________________________________________________
vector<float> XSM_MLP::CreateTMVAInput(IsotopicVector isotopicvector,int TimeStep)
{
    DBGL
    vector<float> tmva_input;

    float Time = 0;

    IsotopicVector IVInputTMVA;
    map<ZAI , string >::iterator it_ZAI_s;

    for ( it_ZAI_s = fMapOfTMVAVariableNames.begin()  ; it_ZAI_s != fMapOfTMVAVariableNames.end() ; it_ZAI_s++)
    {
        IVInputTMVA +=  ((*it_ZAI_s).first) * 1;
    }
    
    // build IV containing only the ZAI known by the Model
    IsotopicVector IVAccordingToUserInfoFile    = isotopicvector.GetThisComposition(IVInputTMVA);
    
    // Normalize the vector
    double Ntot                     = IVAccordingToUserInfoFile.GetSumOfAll();
    IVAccordingToUserInfoFile           = IVAccordingToUserInfoFile / Ntot;

    // Add value in the input vector
    for ( it_ZAI_s = fMapOfTMVAVariableNames.begin() ; it_ZAI_s != fMapOfTMVAVariableNames.end() ; it_ZAI_s++)
    {
        tmva_input.push_back( IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it_ZAI_s).first ) );
    }

    tmva_input.push_back(fMLP_Time[TimeStep]);

    return tmva_input;
}
//________________________________________________________________________
EvolutionData XSM_MLP::GetCrossSections(IsotopicVector IV, double t) {
    DBGL

    string dir = fTMVAWeightFolder;
    if (dir[dir.size() - 1] != '/') {
        dir += "/";
    }

    EvolutionData EvolutionDataFromMLP = EvolutionData();

    map<ZAI, TGraph*> ExtrapolatedXS[3];
    /*************DATA BASE INFO****************/
    EvolutionDataFromMLP.SetReactorType(fDBRType);
    EvolutionDataFromMLP.SetFuelType(fDBFType);
    EvolutionDataFromMLP.SetPower(fDBPower);
    EvolutionDataFromMLP.SetHeavyMetalMass(fDBHMMass);
    /************* The Cross sections***********/
    int Z[int(fWeightFiles.size())];
    int A[int(fWeightFiles.size())];
    int I[int(fWeightFiles.size())];
    int Reaction[int(fWeightFiles.size())];

    for (int i = 0; i < int(fWeightFiles.size()); i++) {
        Z[i] = -2;
        A[i] = -2;
        I[i] = -2;
        Reaction[i] = -2;
        ReadWeightFile(fWeightFiles[i], Z[i], A[i], I[i], Reaction[i]);
    }
    
    vector<float> tmva_input = CreateTMVAInput(IV, 0);
    for (int TimeStep = 0; TimeStep < int(fMLP_Time.size()); TimeStep++) {
      tmva_input.back() = fMLP_Time[TimeStep];  
      for (int i = 0; i < int(fWeightFiles.size()); i++) {
            if (Z[i] >= GetZAIThreshold()) {
                fReader[i]->SetInputData(tmva_input);
                double XSValue = fReader[i]->EvaluateRegression("MLP method")[0];
                // std::cout << "XSValue " << XSValue << std::endl;
                pair<map<ZAI, TGraph*>::iterator, bool> IResult =
                        ExtrapolatedXS[Reaction[i]].insert(
                                pair<ZAI, TGraph*>(ZAI(Z[i], A[i], I[i]), new TGraph()));

                if (IResult.second) {
                    (IResult.first)
                            ->second->SetPoint(0, (double)fMLP_Time[TimeStep], XSValue);

                } else {
                    (IResult.first)
                            ->second->SetPoint((IResult.first)->second->GetN(),
                                                                 (double)fMLP_Time[TimeStep], XSValue);
                }
            }
        }
    }

    /**********Sorting TGraph*********/
    for (int x = 0; x < 3; x++) {
        map<ZAI, TGraph*>::iterator it;
        for (it = ExtrapolatedXS[x].begin(); it != ExtrapolatedXS[x].end(); it++)
            it->second->Sort();
    }
    /**********Filling Matrices*/
    EvolutionDataFromMLP.SetFissionXS(ExtrapolatedXS[0]);
    EvolutionDataFromMLP.SetCaptureXS(ExtrapolatedXS[1]);
    EvolutionDataFromMLP.Setn2nXS(ExtrapolatedXS[2]);

    DBGL return EvolutionDataFromMLP;
}

//________________________________________________________________________
