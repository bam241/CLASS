/*
this code take a fuel composition and print the BU Max
=> This is for MOX EUS Model only
*/
#include "CLASSHeaders.hxx"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <algorithm>
#include <vector>
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
  float U5,U8,Pu8,Pu9,Pu10,Pu11,Pu12,Am1, Time;

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
  FuelToTest.Add(92,235,0,atof(argv[1])); U5   = atof(argv[1]);
  FuelToTest.Add(92,238,0,atof(argv[2])); U8   = atof(argv[2]);
  FuelToTest.Add(94,238,0,atof(argv[3])); Pu8  = atof(argv[3]);
  FuelToTest.Add(94,239,0,atof(argv[4])); Pu9  = atof(argv[4]);
  FuelToTest.Add(94,240,0,atof(argv[5])); Pu10 = atof(argv[5]);
  FuelToTest.Add(94,241,0,atof(argv[6])); Pu11 = atof(argv[6]);
  FuelToTest.Add(94,242,0,atof(argv[7])); Pu12 = atof(argv[7]);
  FuelToTest.Add(95,241,0,atof(argv[8])); Am1  = atof(argv[8]);
  Time=atof(argv[9]) * (cYear);

// Reader
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  reader->AddVariable( "U5"  ,&U5);
  reader->AddVariable( "U8"  ,&U8);
  reader->AddVariable( "Pu8" ,&Pu8);
  reader->AddVariable( "Pu9" ,&Pu9);
  reader->AddVariable( "Pu10",&Pu10);
  reader->AddVariable( "Pu11",&Pu11);
  reader->AddVariable( "Pu12",&Pu12);
  reader->AddVariable( "Am1" ,&Am1);
  reader->AddVariable( "Time",&Time);
  reader->BookMVA("MLP method", XML_E);

// Kinf
  double Kinf_t = (reader->EvaluateRegression("MLP method"))[0]; cout<<Kinf_t<<endl;
//BU
  double BU = E_EUS->CalculateBurnUpMax(FuelToTest, ModelParameter);

}
/*
g++ -std=c++11 -o PWR_MOX_PrintMLP PWR_MOX_PrintMLP.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -lTMVA

PWR_MOX_PrintMLP 0.00237622 0.948111 0.00346011 0.0224409 0.0101809 0.0075659 0.00475449 0.0011101 0
PWR_MOX_PrintMLP 0.00230808 0.920924 0.00461033 0.0277897 0.0184003 0.0108105 0.012515 0.00264214 0
*/
