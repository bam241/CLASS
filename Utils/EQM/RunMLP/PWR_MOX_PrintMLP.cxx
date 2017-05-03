/*
this code take a fuel composition and print the BU Max
=> This is for MOX EUS Model only
*/
#include "CLASSReader.hxx"
#include "CLASSHeaders.hxx"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
/*
#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "IsotopicVector.hxx"
#include <math.h>
#include "TTree.h"
#include <map>
#include "CLASSObject.hxx"
*/
using namespace std;

//---------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------- MAIN ------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  CLASSLogger *Logger = new CLASSLogger("C_LOG.log",0,2);

// Parameters
  double k_th_E = 1.034;
  int N_batches = 3;
  double T; //time

  // Path
  string CLASS_PATH   = getenv("CLASS_PATH");
  string PATH_TO_DATA = CLASS_PATH + "/DATA_BASES/";
  string NFO_E        = PATH_TO_DATA + "PWR/MOX_EUS/EQModel/NFO/PWR_MOxEUS_Pu_0_16.nfo";
  string XML_E        = PATH_TO_DATA + "PWR/MOX_EUS/EQModel/XML/PWR_MOxEUS_Pu_0_16.xml";


  // EQM Model
  EquivalenceModel* E_EUS = new EquivalenceModel(Logger, XML_E,  NFO_E );
  E_EUS->SetModelParameter("kThreshold",  k_th_E);
  E_EUS->SetModelParameter("NumberOfBatch",  N_batches);

  // Model parameter
  map<string, double> ModelParameter;
  ModelParameter["kThreshold"]  = k_th_E;
  ModelParameter["NumberOfBatch"] = N_batches;

  // Fuel To Test
  IsotopicVector FuelToTest;
  FuelToTest.Add(92,235,0,atof(argv[1])); //U5
  FuelToTest.Add(92,238,0,atof(argv[2])); //U8
  FuelToTest.Add(94,238,0,atof(argv[3])); //Pu8
  FuelToTest.Add(94,239,0,atof(argv[4])); //Pu9
  FuelToTest.Add(94,240,0,atof(argv[5])); //Pu10
  FuelToTest.Add(94,241,0,atof(argv[6])); //Pu11
  FuelToTest.Add(94,242,0,atof(argv[7])); //Pu12
  FuelToTest.Add(95,241,0,atof(argv[8])); //Am1

  map<ZAI,string> MapOfTMVAVariableNames;
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(92,235,0), "U5"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(92,238,0), "U8"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(94,238,0), "Pu8"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(94,239,0), "Pu9"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(94,240,0), "Pu10"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(94,241,0), "Pu11"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(94,242,0), "Pu12"));
  MapOfTMVAVariableNames.insert(pair<ZAI,string>(ZAI(95,241,0), "Am1"));
/*
  // Reader
  
  CLASSReader * reader = new CLASSReader(MapOfTMVAVariableNames);
  reader->AddVariable("Time");
  reader->BookMVA("MLP method", XML_E);
  TTree* InputTree = E_EUS->CreateTMVAInputTree(FuelToTest,T);
  reader->SetInputData(InputTree);
  
// Kinf
  double Kinf_t = reader->EvaluateRegression( "MLP method" )[0];
*/
//BU
  double BU = E_EUS->CalculateBurnUpMax(FuelToTest, ModelParameter);


//  delete InputTree;
}
/*
g++ -std=c++11 -o PWR_MOX_PrintMLP PWR_MOX_PrintMLP.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp 
*/

 // 0.00237622 0.948111 0.00346011 0.0224409 0.0101809 0.0075659 0.00475449 0.0011101
 // 0.00230808 0.920924 0.00461033 0.0277897 0.0184003 0.0108105 0.012515 0.00264214
