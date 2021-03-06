#include "EQ_OneParameter.hxx"
#include "CLASSMethod.hxx"
#include "StringLine.hxx"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TSystem.h"

//________________________________________________________________________
//________________________________________________________________________
EQ_OneParameter::EQ_OneParameter(string TMVAXMLFilePath, string TMVANFOFilePath)
    : EquivalenceModel(new CLASSLogger("EQ_OneParameter.log")) {
  fUseTMVAPredictor = true;

  fMaxIterration = 500;
  fPCMprecision = 10;

  fTMVAXMLFilePath = TMVAXMLFilePath;
  fTMVANFOFilePath = TMVANFOFilePath;

  fDBRType = "";
  fDBFType = "";
  fSpecificPower = 0;
  fMaximalBU = 0;
  fTargetParameter = "";
  fTargetParameterStDev = 0;
  fBuffer = "";
  fPredictorType = "";
  fOutput = "";

  LoadKeyword();  // Load Key words defineds in NFO file
  ReadNFO();      // Getting information from file NFO

  BookTMVAReader();
  CheckParameterValidity();

  INFO << "__An equivalence model has been define__" << endl;
  INFO << "\tThe TMVA weights file is :" << fTMVAXMLFilePath << endl;
  INFO << "\tThe TMVA NFO file is :" << fTMVANFOFilePath << endl;
  PrintInfo();
}
//________________________________________________________________________
EQ_OneParameter::EQ_OneParameter(CLASSLogger *log, string TMVAXMLFilePath,
                                 string TMVANFOFilePath)
    : EquivalenceModel(log) {
  fUseTMVAPredictor = true;

  fMaxIterration = 500;
  fread = false;
  fPCMprecision = 10;

  fTMVAXMLFilePath = TMVAXMLFilePath;
  fTMVANFOFilePath = TMVANFOFilePath;

  fDBRType = "";
  fDBFType = "";
  fSpecificPower = 0;
  fMaximalBU = 0;
  fTargetParameter = "";
  fTargetParameterStDev = 0;
  fBuffer = "";
  fPredictorType = "";
  fOutput = "";

  LoadKeyword();  // Load Key words defineds in NFO file
  ReadNFO();      // Getting information from file NFO

  BookTMVAReader();
  CheckParameterValidity();

  INFO << "__An equivalence model has been define__" << endl;
  INFO << "\tThe TMVA weights file is :" << fTMVAXMLFilePath << endl;
  INFO << "\tThe TMVA NFO file is :" << fTMVANFOFilePath << endl;
  PrintInfo();
}
//________________________________________________________________________
EQ_OneParameter::~EQ_OneParameter() { delete fReader; }

void EQ_OneParameter::CheckParameterValidity() {
  // Check if any information is missing in NFO file
  bool NFOValidity = true;

  if (fZAILimits.empty()) {
    ERROR << "Missing information for k_zail in : " << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (fDBRType.empty()) {
    ERROR << "Missing information for k_reactor in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (fDBFType.empty()) {
    ERROR << "Missing information for k_fuel in : " << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (!fSpecificPower) {
    ERROR << "Missing information for k_specpower in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (!fMaximalBU) {
    ERROR << "Missing information for k_maxburnup in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (fStreamListEqMMassFractionMin.empty() ||
      fStreamListEqMMassFractionMax.empty()) {
    ERROR << "Missing information for k_massfractionmin and/or "
             "k_massfractionmax in : "
          << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (fStreamList.empty()) {
    ERROR << "Missing information for k_list in : " << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (fMapOfTMVAVariableNames.empty()) {
    ERROR << "Missing information for k_zainame in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (fTargetParameter.empty()) {
    ERROR << "Missing information for k_targetparameter in : "
          << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (!fTargetParameterStDev) {
    ERROR << "Missing information for fTargetParameterStDev in : "
          << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (fModelParameter.empty()) {
    ERROR << "Missing information for k_modelparameter in : "
          << fTMVANFOFilePath << endl;
    NFOValidity = false;
  }
  if (fBuffer.empty()) {
    ERROR << "Missing information for k_buffer in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (fPredictorType.empty()) {
    ERROR << "Missing information for k_predictortype in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (fOutput.empty()) {
    ERROR << "Missing information for k_output in : " << fTMVANFOFilePath
          << endl;
    NFOValidity = false;
  }
  if (!NFOValidity) {
    exit(1);
  }
}
void EQ_OneParameter::BookTMVAReader() {
  DBGL;
  fReader = new CLASSReader();

  for (int i = 0; i < fListOfNonZaiTMVAVariables.size(); i++) {
    fReader->AddVariable(fListOfNonZaiTMVAVariables[i].second);
  }

  map<ZAI, string>::iterator it;
  for (it = fMapOfTMVAVariableNames.begin();
       it != fMapOfTMVAVariableNames.end(); it++) {
    fReader->AddVariable(it->second);
  }

  if (fTargetParameter == "BurnUpMax" || fDBRType == "SFR") {
    fReader->AddVariable("Time");
  }
  fReader->BookMVA("MLP method", fTMVAXMLFilePath);
  DBGL;
}

//________________________________________________________________________
IsotopicVector EQ_OneParameter::BuildFuelToTest(
    map<string, vector<double> > &lambda,
    map<string, vector<IsotopicVector> > const &StreamArray, double HMMass,
    map<string, bool> StreamListIsBuffer) {
  DBGL;
  // Iterators declaration
  map<string, vector<IsotopicVector> >::const_iterator it_s_vIV;
  map<string, bool>::iterator it_s_B;

  // Find the buffer and set its lambda to 0
  string BufferDenomination = "";
  for (it_s_B = StreamListIsBuffer.begin(); it_s_B != StreamListIsBuffer.end();
       it_s_B++) {
    if ((*it_s_B).second == true) {
      BufferDenomination = (*it_s_B).first;
    }
  }

  for (int i = 0; i < lambda[BufferDenomination].size(); i++) {
    lambda[BufferDenomination][i] = 0;
  }

  // Build an IV with all materials besides buffer to get the total mass of
  // others materials
  IsotopicVector IV;
  for (it_s_vIV = StreamArray.begin(); it_s_vIV != StreamArray.end();
       it_s_vIV++) {
    for (int i = 0; i < (int)StreamArray.at(it_s_vIV->first).size(); i++) {
      IV += lambda[(*it_s_vIV).first][i] * StreamArray.at(it_s_vIV->first)[i];
    }
  }

  // Calculate MassBuffer
  double MassBuffer = HMMass - IV.GetTotalMass() * 1e06;

  // Set buffer lambda according to MassBuffer
  ConvertMassToLambdaVector(BufferDenomination, lambda[BufferDenomination],
                            MassBuffer, StreamArray.at(BufferDenomination));

  IV.Clear();

  // Build fuel with all materials
  for (it_s_vIV = StreamArray.begin(); it_s_vIV != StreamArray.end();
       it_s_vIV++) {
    for (int i = 0; i < (int)StreamArray.at(it_s_vIV->first).size(); i++) {
      IV += lambda[(*it_s_vIV).first][i] * StreamArray.at(it_s_vIV->first)[i];
    }
  }
  DBGL;
  return IV;
}

//________________________________________________________________________
map<string, vector<double> > EQ_OneParameter::BuildFuel(
    double BurnUp, double HMMass,
    map<string, vector<IsotopicVector> > StreamArray,
    map<string, double> StreamListFPMassFractionMin,
    map<string, double> StreamListFPMassFractionMax,
    map<int, string> StreamListPriority, map<string, bool> StreamListIsBuffer) {
  DBGL;

  HMMass *= 1e6;  // Unit conversion : tons to gram

  map<string, vector<double> > lambda;  // map containing name of the list and
                                        // associated vector of proportions
                                        // taken from stocks
  // Iterators declaration
  map<string, vector<IsotopicVector> >::iterator it_s_vIV;
  map<string, vector<double> >::iterator it_s_vD;
  map<string, IsotopicVector>::iterator it_s_IV;
  map<string, double>::iterator it_s_D;
  map<int, string>::iterator it_i_s;

  // Initialize lambda to 0 //
  for (it_s_vIV = StreamArray.begin(); it_s_vIV != StreamArray.end();
       it_s_vIV++) {
    for (int i = 0; i < (int)StreamArray[(*it_s_vIV).first].size(); i++) {
      lambda[(*it_s_vIV).first].push_back(0);
    }
  }

  // Test if there is at least one stock available in each list, otherwise fuel
  // is not built //
  bool BreakReturnLambda = true;
  for (it_s_vIV = StreamArray.begin(); it_s_vIV != StreamArray.end();
       it_s_vIV++) {
    if (StreamArray[(*it_s_vIV).first].size() != 0) {
      BreakReturnLambda &= false;
    }
  }

  if (BreakReturnLambda) {
    WARNING << " No stock available for stream : " << (*it_s_vIV).first
            << ".  Fuel not built." << endl;
    SetLambdaToErrorCode(lambda[(*it_s_vIV).first]);
    return lambda;
  }

  // Check if the targeted burn-up is lower than maximum burn-up of model //
  if (BurnUp > fMaximalBU) {
    ERROR << " Targeted burn-up is higher than maximum burn-up defined in NFO "
             "file..."
          << endl;
    ERROR << " Targeted burn-up : " << BurnUp << " GWd/t" << endl;
    ERROR << " Maximum burn-up : " << fMaximalBU << " GWd/t" << endl;
    exit(1);
  }

  // Fissile fraction calculation is needed.
  if (fUseTMVAPredictor) {
    // Check if EQ_OneParameter->SetTMVAXMLFilePath() and/or
    // EQ_OneParameter->SetTMVANFOFilePath() have been defined
    if (fTMVAXMLFilePath.empty() || fTMVANFOFilePath.empty()) {
      ERROR << " TMVA XML and/or NFO File path are not defined..." << endl;
      ERROR << " You have to use EQ_OneParameter->SetTMVAXMLFilePath() and/or "
               "EQ_OneParameter->SetTMVANFOFilePath() methods."
            << endl;
      exit(1);
    }

    double TargetParameterValue = 0;

    for (it_s_D = fModelParameter.begin(); it_s_D != fModelParameter.end();
         it_s_D++) {
      if (fModelParameter[(*it_s_D).first] == -1) {
        ERROR << "Model parameter ( " << fModelParameter[(*it_s_D).first]
              << " ) value is not defined in the input." << endl;
        ERROR << "Use EqM->SetModelParameter( \" " << (*it_s_D).first
              << " \", value) to define it." << endl;
        exit(1);
      }
    }

    if (fTargetParameter == "BurnUpMax") {
      TargetParameterValue = BurnUp;
    } else if (fTargetParameter == "keffBOC") {
      TargetParameterValue = fModelParameter["keffBOC"];
    } else {
      ERROR << "Target parameter defined in InformationFile ( "
            << fTargetParameter << " ) doesn't exist." << endl;
      ERROR << "Possible target parameters for the moment are : " << endl;
      ERROR << " - BurnUpMax - Used for PWR" << endl;
      ERROR << " - keffBOC - Used for SFR" << endl;
      exit(1);
    }

    /// Search for the minimum and maximum fraction of each material in fuel ///
    map<string, double> StreamListMassFractionMin;
    map<string, double> StreamListMassFractionMax;
    for (it_s_D = StreamListFPMassFractionMin.begin();
         it_s_D != StreamListFPMassFractionMin.end(); it_s_D++) {
      if (StreamListFPMassFractionMin[(*it_s_D).first] <
          fStreamListEqMMassFractionMin[(*it_s_D).first])  // if limits FP are
                                                           // lower than limits
                                                           // EqM
      {
        ERROR << " User mass fraction min requirement is lower than the model "
                 "mass fraction min for list  : "
              << (*it_s_D).first << endl;
        ERROR << " User mass fraction min requirement : "
              << StreamListFPMassFractionMin[(*it_s_D).first] << endl;
        ERROR << " Model mass fraction min requirement : "
              << fStreamListEqMMassFractionMin[(*it_s_D).first] << endl;
        exit(1);
      } else {
        StreamListMassFractionMin[(*it_s_D).first] =
            StreamListFPMassFractionMin[(*it_s_D).first];
      }
    }

    for (it_s_D = StreamListFPMassFractionMax.begin();
         it_s_D != StreamListFPMassFractionMax.end(); it_s_D++) {
      if (StreamListFPMassFractionMax[(*it_s_D).first] >
          fStreamListEqMMassFractionMax[(*it_s_D).first])  // if limits FP are
                                                           // higher than limits
                                                           // EqM
      {
        ERROR << " User mass fraction max requirement is higher than the model "
                 "mass fraction max for list  : "
              << (*it_s_D).first << endl;
        ERROR << " User mass fraction max requirement : "
              << StreamListFPMassFractionMax[(*it_s_D).first] << endl;
        ERROR << " Model mass fraction max requirement : "
              << fStreamListEqMMassFractionMax[(*it_s_D).first] << endl;
        exit(1);
      } else {
        StreamListMassFractionMax[(*it_s_D).first] =
            StreamListFPMassFractionMax[(*it_s_D).first];
      }
    }

    // Calculate Total mass in stock for each stream and fill fTotalMassInStocks
    StocksTotalMassCalculation(StreamArray);

    // Check if there is enough material in stock to satisfy mass fraction min
    // //
    BreakReturnLambda = false;
    for (it_s_D = StreamListMassFractionMin.begin();
         it_s_D != StreamListMassFractionMin.end(); it_s_D++) {
      if (fTotalMassInStocks[(*it_s_D).first] <
          HMMass * StreamListMassFractionMin[(*it_s_D).first]) {
        WARNING << " Not enough material  : " << (*it_s_D).first
                << " in stocks to reach the build fuel lower limit of "
                << StreamListMassFractionMin[(*it_s_D).first]
                << " reactor mass.  Fuel not built." << endl;
        SetLambdaToErrorCode(lambda[(*it_s_D).first]);
        BreakReturnLambda = true;
      }
    }
    if (BreakReturnLambda) {
      return lambda;
    }

    // Check if there is enough material in stock to satisfy mass fraction max,
    // if not mass fraction max is set to MassINStock/MassReactor//
    for (it_s_D = StreamListMassFractionMax.begin();
         it_s_D != StreamListMassFractionMax.end(); it_s_D++) {
      if (fTotalMassInStocks[(*it_s_D).first] <
          HMMass * StreamListMassFractionMax[(*it_s_D).first]) {
        StreamListMassFractionMax[(*it_s_D).first] =
            fTotalMassInStocks[(*it_s_D).first] / HMMass;
        WARNING << " Not enough material  : " << (*it_s_D).first
                << " in stocks to reach the build fuel higher limit of "
                << StreamListMassFractionMax[(*it_s_D).first]
                << " reactor mass. " << endl;
        WARNING << " Mass fraction max of material :  " << (*it_s_D).first
                << " is set to MassInStock/HMMassReactor : "
                << StreamListMassFractionMax[(*it_s_D).first] << endl;
      }
    }

    // Check if TargetParameter is inside [TargetParameterMin,
    // TargetParameterMax] associated to fraction Min et Max//

    map<string, double> MassMin;
    map<string, double> MassMax;

    map<string, double> TargetParameterMin;
    map<string, double> TargetParameterMax;

    IsotopicVector FuelToTest;

    bool TargetParameterIncluded = false;
    for (it_i_s = StreamListPriority.begin();
         it_i_s != StreamListPriority.end(); it_i_s++) {
      // Calculate TargetParameterMin for each possibility : min1 ; max1 + min2
      // ;  max1 + max2 + min3 ....
      MassMin[(*it_i_s).second] =
          HMMass * StreamListMassFractionMin[(*it_i_s).second];
      ConvertMassToLambdaVector((*it_i_s).second, lambda[(*it_i_s).second],
                                MassMin[(*it_i_s).second],
                                StreamArray[(*it_i_s).second]);
      FuelToTest =
          BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
      FuelToTest = FuelToTest / FuelToTest.GetSumOfAll();
      TargetParameterMin[(*it_i_s).second] =
          CalculateTargetParameter(FuelToTest, fTargetParameter);

      // Check is TargetParameterMin < TargetParameter
      if (TargetParameterMin[(*it_i_s).second] > TargetParameterValue) {
        if ((*it_i_s).first == 1)  // Minimum of first material is too high
        {
          WARNING << "CRITICAL ! Minimum parameter value associated to the "
                     "first priority material ( "
                  << (*it_i_s).second << " ) is higher than targeted parameter."
                  << endl;
          WARNING << "Targeted parameter : " << fTargetParameter << " = "
                  << TargetParameterValue << endl;
          WARNING << "Minimum parameter value : "
                  << TargetParameterMin[(*it_i_s).second] << endl;
          WARNING << "Try to increase targeted parameter." << endl;
          SetLambdaToErrorCode(lambda[(*it_i_s).second]);
          return lambda;
          DBGL;
        } else if ((*it_i_s).first >
                   1)  // TargetParameter is located between max n-1 and min n
        {
          WARNING << "CRITICAL ! Targeted parameter value ( "
                  << fTargetParameter << " ) is located between 2 materials. "
                  << endl;
          it_i_s--;
          WARNING << fTargetParameter
                  << " of max fraction of material : " << (*it_i_s).second
                  << " ---> " << TargetParameterMax[(*it_i_s).second] << endl;
          it_i_s++;
          WARNING << fTargetParameter
                  << " of min fraction of material : " << (*it_i_s).second
                  << " ---> " << TargetParameterMin[(*it_i_s).second] << endl;
          WARNING << "Targeted " << fTargetParameter << " : "
                  << TargetParameterValue << endl;
          WARNING << "Try to decrease mimimum fraction of : "
                  << (*it_i_s).second << endl;
          SetLambdaToErrorCode(lambda[(*it_i_s).second]);
          return lambda;
        }
      }
      FuelToTest.Clear();

      // Calculate TargetParameter max for each possibility : max1 ; max1 + max2
      // ;  max1 + max2 + max3 ....
      MassMax[(*it_i_s).second] =
          HMMass * StreamListMassFractionMax[(*it_i_s).second];
      ConvertMassToLambdaVector((*it_i_s).second, lambda[(*it_i_s).second],
                                MassMax[(*it_i_s).second],
                                StreamArray[(*it_i_s).second]);
      FuelToTest =
          BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
      FuelToTest = FuelToTest / FuelToTest.GetSumOfAll();
      TargetParameterMax[(*it_i_s).second] =
          CalculateTargetParameter(FuelToTest, fTargetParameter);

      if (TargetParameterMax[(*it_i_s).second] >= TargetParameterValue) {
        TargetParameterIncluded = true;
        break;
      }
    }

    // Check if target parameter increases monotously with the material mass
    CheckTargetParameterConsistency(StreamListPriority, TargetParameterMin,
                                    TargetParameterMax);

    if (!TargetParameterIncluded) {
      WARNING << "CRITICAL ! Maximum reachable " << fTargetParameter
              << " is lower than targeted " << fTargetParameter << ". " << endl;
      WARNING << "Targeted " << fTargetParameter << " = "
              << TargetParameterValue << endl;
      WARNING << "Maximum reachable " << fTargetParameter << " : "
              << TargetParameterMax[(*--StreamListPriority.end()).second]
              << endl;
      WARNING << "Try to increase maximum fraction of materials, or decrease "
              << fTargetParameter << " ." << endl;
      SetLambdaToErrorCode(lambda[(*--StreamListPriority.end()).second]);
      return lambda;
    }

    // Search the TargetParameterValue location in the mass damain //
    string MaterialToSearch = (*it_i_s).second;
    double CalculatedTargetParameter =
        TargetParameterMax[MaterialToSearch];  // Algo start with maximum point
    double MassToAdd =
        MassMax[MaterialToSearch];  // Algo start with maximum point

    double LastMassMinus =
        MassMin[MaterialToSearch];  // Used in bissection method
    double LastMassPlus =
        MassMax[MaterialToSearch];  // Used in bissection method

    int count = 0;

    FuelToTest.Clear();

    do {
      if (count > fMaxIterration) {
        ERROR << "CRITICAL ! Can't manage to predict fissile content\nHint : "
                 "Try to decrease the precision on the target parameter using "
                 ":\nYourEQ_OneParameter->SetTargetParameterStDev(Precision); "
              << endl;
        ERROR << "Targeted " << fTargetParameter << " : "
              << TargetParameterValue << endl;
        ERROR << "Last calculated " << fTargetParameter << " : "
              << CalculatedTargetParameter << endl;
        ERROR << "Last Fresh fuel normalized composition : " << endl;
        ERROR << FuelToTest.sPrint() << endl;
        exit(1);
      }

      if ((CalculatedTargetParameter - TargetParameterValue) <
          0)  // Need to add more fissile material in fuel
      {
        LastMassMinus = MassToAdd;
        MassToAdd = MassToAdd + fabs(LastMassPlus - MassToAdd) / 2.;
      } else if ((CalculatedTargetParameter - TargetParameterValue) >
                 0)  // Need to add less fissile material in fuel
      {
        LastMassPlus = MassToAdd;
        MassToAdd = MassToAdd - fabs(LastMassMinus - MassToAdd) / 2.;
      }
      ConvertMassToLambdaVector(MaterialToSearch, lambda[MaterialToSearch],
                                MassToAdd, StreamArray[MaterialToSearch]);
      FuelToTest =
          BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
      FuelToTest = FuelToTest / FuelToTest.GetSumOfAll();
      CalculatedTargetParameter =
          CalculateTargetParameter(FuelToTest, fTargetParameter);

      count++;

    } while (fabs(TargetParameterValue - CalculatedTargetParameter) >
             GetTargetParameterStDev() * TargetParameterValue);
  }
  // Fissile fraction is imposed by the FP
  // No need to use algo
  // Simplified fuel building
  else {
    // Check if EQ_OneParameter->SetTMVANFOFilePath() have been defined
    if (fTMVANFOFilePath.empty()) {
      ERROR << " TMVA NFO File path is not defined..." << endl;
      ERROR << " You have to use EQ_OneParameter->SetTMVANFOFilePath() methods."
            << endl;
      exit(1);
    }

    /// Search for the  fraction of each material in fuel ///
    map<string, double> StreamListMassFraction;
    for (it_s_D = StreamListFPMassFractionMin.begin();
         it_s_D != StreamListFPMassFractionMin.end(); it_s_D++) {
      if (StreamListFPMassFractionMin[(*it_s_D).first] <
          fStreamListEqMMassFractionMin[(*it_s_D).first])  // if limits FP are
                                                           // lower than limits
                                                           // EqM
      {
        ERROR << " User mass fraction requirement is lower than the model mass "
                 "fraction min for list  : "
              << (*it_s_D).first << endl;
        ERROR << " User mass fraction requirement : "
              << StreamListFPMassFractionMin[(*it_s_D).first] << endl;
        ERROR << " Model mass fraction min requirement : "
              << fStreamListEqMMassFractionMin[(*it_s_D).first] << endl;
        exit(1);
      } else if (StreamListFPMassFractionMax[(*it_s_D).first] >
                 fStreamListEqMMassFractionMax[(*it_s_D).first])  // if limits
                                                                  // FP are
                                                                  // higher than
                                                                  // limits EqM
      {
        ERROR << " User mass fraction requirement is higher than the model "
                 "mass fraction max for list  : "
              << (*it_s_D).first << endl;
        ERROR << " User mass fraction requirement : "
              << StreamListFPMassFractionMax[(*it_s_D).first] << endl;
        ERROR << " Model mass fraction max requirement : "
              << fStreamListEqMMassFractionMax[(*it_s_D).first] << endl;
        exit(1);
      } else {
        StreamListMassFraction[(*it_s_D).first] =
            StreamListFPMassFractionMin[(*it_s_D).first];  // Because here, min
                                                           // = max
      }
    }

    // Calculate Total mass in stock for each stream and fill fTotalMassInStocks
    StocksTotalMassCalculation(StreamArray);

    // Check if there is enough material in stock to satisfy requested mass
    // fraction //
    BreakReturnLambda = false;
    for (it_s_D = StreamListMassFraction.begin();
         it_s_D != StreamListMassFraction.end(); it_s_D++) {
      if (fTotalMassInStocks[(*it_s_D).first] <
          HMMass * StreamListMassFraction[(*it_s_D).first]) {
        WARNING << " Not enough material  : " << (*it_s_D).first
                << " in stocks to reach the build fuel limit of "
                << StreamListMassFraction[(*it_s_D).first]
                << " reactor mass.  Fuel not built." << endl;
        SetLambdaToErrorCode(lambda[(*it_s_D).first]);
        BreakReturnLambda = true;
      }
    }
    if (BreakReturnLambda) {
      return lambda;
    }

    IsotopicVector FuelToTest;
    // Build Fuel
    for (it_i_s = StreamListPriority.begin();
         it_i_s != StreamListPriority.end(); it_i_s++) {
      FuelToTest.Clear();
      ConvertMassToLambdaVector(
          (*it_i_s).second, lambda[(*it_i_s).second],
          HMMass * StreamListMassFraction[(*it_i_s).second],
          StreamArray[(*it_i_s).second]);
      FuelToTest =
          BuildFuelToTest(lambda, StreamArray, HMMass, StreamListIsBuffer);
    }
  }

  // Final builded fuel
  IsotopicVector IVStream;
  for (it_s_vD = lambda.begin(); it_s_vD != lambda.end(); it_s_vD++) {
    for (int i = 0; i < (int)lambda[(*it_s_vD).first].size(); i++) {
      IVStream +=
          lambda[(*it_s_vD).first][i] * StreamArray[(*it_s_vD).first][i];
    }
  }

  // Check if BuildedFuel is in Model isotopic bounds
  (*this).isIVInDomain(IVStream);

  for (it_s_vD = lambda.begin(); it_s_vD != lambda.end(); it_s_vD++) {
    DBGV("Lambda vector : " << (*it_s_vD).first);

    for (int i = 0; i < (int)lambda[(*it_s_vD).first].size(); i++) {
      DBGV(lambda[(*it_s_vD).first][i]);
    }
  }

  DBGL;
  return lambda;
}

//________________________________________________________________________
vector<float> EQ_OneParameter::CreateTMVAInput(IsotopicVector TheFreshfuel,
                                               double ThisTime) {
  /******Create Input data tree to be interpreted by TMVA::Reader***/
  vector<float> tmva_input;

  float Time = 0;

  IsotopicVector IVInputTMVA;
  map<ZAI, string>::iterator it_ZAI_s;

  for (int j = 0; j < fListOfNonZaiTMVAVariables.size(); j++) {
    tmva_input.push_back(fListOfNonZaiTMVAVariables[j].first);
  }

  for (it_ZAI_s = fMapOfTMVAVariableNames.begin();
       it_ZAI_s != fMapOfTMVAVariableNames.end(); it_ZAI_s++) {
    IVInputTMVA += ((*it_ZAI_s).first) * 1;
  }

  // build IV containing only the ZAI known by the Model
  IsotopicVector IVAccordingToUserInfoFile =
      TheFreshfuel.GetThisComposition(IVInputTMVA);

  // Normalize the vector
  double Ntot = IVAccordingToUserInfoFile.GetSumOfAll();
  IVAccordingToUserInfoFile = IVAccordingToUserInfoFile / Ntot;

  // Add value in the input vector
  for (it_ZAI_s = fMapOfTMVAVariableNames.begin();
       it_ZAI_s != fMapOfTMVAVariableNames.end(); it_ZAI_s++) {
    tmva_input.push_back(
        IVAccordingToUserInfoFile.GetZAIIsotopicQuantity((*it_ZAI_s).first));
  }

  if (ThisTime != -1) tmva_input.push_back(ThisTime);

  return tmva_input;
}
//________________________________________________________________________
void EQ_OneParameter::CheckTargetParameterConsistency(
    map<int, string> StreamListPriority, map<string, double> TargetParameterMin,
    map<string, double> TargetParameterMax) {
  map<int, string>::iterator it_i_s;

  // Loop on priority order to check if target parameter increases monotously
  // with the material mass
  for (it_i_s = StreamListPriority.begin(); it_i_s != StreamListPriority.end();
       it_i_s++) {
    double TargetParameterUp =
        -1.0;  // to be sure BUMin is > to BUmax even if BUmin is zero
    double TargetParameterDown = 0.0;

    if (TargetParameterMin.find((*it_i_s).second) == TargetParameterMin.end()) {
      break;  // if material is not in map, break the loop
    }
    TargetParameterDown = TargetParameterMin[(*it_i_s).second];

    if (TargetParameterDown < 0.0) {
      ERROR << "Target parameter evolution should always be positive." << endl;
      ERROR << "TargetParameterDown = " << TargetParameterDown
            << " is negative " << endl;
      ERROR << "Check the evolution..." << endl;
      exit(1);
    }
    TargetParameterUp = TargetParameterMax[(*it_i_s).second];

    if (TargetParameterDown > TargetParameterUp) {
      ERROR << "Target parameter evolution as a function of material mass is "
               "not monotonous."
            << endl;
      ERROR << "TargetParameterDown = " << TargetParameterDown
            << " is greater than  TargetParameterUp = " << TargetParameterUp
            << endl;
      ERROR << "Check the evolution..." << endl;
      exit(1);
    }
  }
}
//________________________________________________________________________
double EQ_OneParameter::CalculateTargetParameter(IsotopicVector TheFuel,
                                                 string TargetParameterName) {
  double ParameterToCalculate = 0;
  if (TargetParameterName == "BurnUpMax")
    ParameterToCalculate = CalculateBurnUpMax(TheFuel);
  else if (TargetParameterName == "keffBOC")
    ParameterToCalculate = CalculateKeffAtBOC(TheFuel);
  else {
    ERROR << "Target parameter defined in InformationFile ( "
          << TargetParameterName << " ) doesn't exist" << endl;
    ERROR << "Possible target parameters for the moment are : BurnUpMax and "
             "keffBOC."
          << endl;
    exit(1);
  }

  return ParameterToCalculate;
}
//________________________________________________________________________
double EQ_OneParameter::CalculateBurnUpMax(IsotopicVector TheFuel) {
  /**************************************************************************/
  // With a dichotomy, the maximal irradiation time (TheFinalTime) is calculated
  // When average Kinf is very close (according "Precision") to the threshold
  // then the corresponding irradiation time is convert in burnup and returned
  /**************************************************************************/
  // Algorithm initialization
  double KThreshold = fModelParameter["kThreshold"];
  int NumberOfBatch = (int)fModelParameter["NumberOfBatch"];
  double OldFinalTimeMinus = 0;
  double MaximumBU = fMaximalBU;
  double MinimumBU = 0;
  double TheFinalTime = BurnupToSecond((MaximumBU - MinimumBU) / 2.);
  double OldFinalTimePlus = BurnupToSecond(MaximumBU);
  double k_av = 0;  // average kinf
  double OldPredictedk_av = 0;

  for (int b = 0; b < NumberOfBatch; b++) {
    float TheTime = (b + 1) * TheFinalTime / NumberOfBatch;

    vector<float> tmva_input = CreateTMVAInput(TheFuel, TheTime);
    fReader->SetInputData(tmva_input);

    OldPredictedk_av += fReader->EvaluateRegression("MLP method")[0];
  }
  OldPredictedk_av /= NumberOfBatch;

  // Algorithm control
  int count = 0;
  int MaximumLoopCount = 500;
  do {
    if (count > MaximumLoopCount) {
      ERROR << "CRITICAL ! Can't manage to predict burnup\nHint : Try to "
               "increase the precision on k effective using :\n "
               "YourEQM_MLP_Kinf->SetPCMPrecision(pcm); with pcm the precision "
               "in pcm (default 10) REDUCE IT\n If this message still appear "
               "mail to leniau@subatech.in2p3.fr\nor "
               "nicolas.thiolliere@subatech.in2p3.fr "
            << endl;
      exit(1);
    }

    if ((OldPredictedk_av - KThreshold) > 0)  // The burnup can be increased
    {
      OldFinalTimeMinus = TheFinalTime;
      TheFinalTime = TheFinalTime + fabs(OldFinalTimePlus - TheFinalTime) / 2.;

      if (SecondToBurnup(TheFinalTime) >=
          (MaximumBU - MaximumBU * GetTargetParameterStDev())) {
        return MaximumBU;
      }
    } else if ((OldPredictedk_av - KThreshold) < 0)  // The burnup is too high
    {
      OldFinalTimePlus = TheFinalTime;
      TheFinalTime = TheFinalTime - fabs(OldFinalTimeMinus - TheFinalTime) / 2.;
      if (SecondToBurnup(TheFinalTime) <
          (MaximumBU - MinimumBU) / 2. * GetTargetParameterStDev()) {
        return 0;
      }
    }

    k_av = 0;
    for (int b = 0; b < NumberOfBatch; b++) {
      float TheTime = (b + 1) * TheFinalTime / NumberOfBatch;
      vector<float> tmva_input = CreateTMVAInput(TheFuel, TheTime);
      fReader->SetInputData(tmva_input);

      k_av += fReader->EvaluateRegression("MLP method")[0];
    }
    k_av /= NumberOfBatch;
    OldPredictedk_av = k_av;
    count++;
  } while (fabs(OldPredictedk_av - KThreshold) > GetPCMPrecision());

  return SecondToBurnup(TheFinalTime);
}

//________________________________________________________________________
double EQ_OneParameter::CalculateKeffAtBOC(IsotopicVector FreshFuel) {
  vector<float> tmva_input = CreateTMVAInput(FreshFuel, -1);
  fReader->SetInputData(tmva_input);

  double keff = fReader->EvaluateRegression("MLP method")[0];

  return keff;
}
//________________________________________________________________________
void EQ_OneParameter::ReadLine(string line) {
  DBGL;
  int pos = 0;
  string keyword = tlc(StringLine::NextWord(line, pos, ' '));

  map<string, EQ_OP_MthPtr>::iterator it = fKeyword.find(keyword);
  if (it != fKeyword.end()) (this->*(it->second))(line);
  DBGL;
}
//________________________________________________________________________
void EQ_OneParameter::LoadKeyword() {
  DBGL;
  EquivalenceModel::LoadKeyword();
  fKeyword.insert(pair<string, EQ_OP_MthPtr>(
      "k_targetparameter", &EQ_OneParameter::ReadTargetParameter));
  fKeyword.insert(pair<string, EQ_OP_MthPtr>(
      "k_modelparameter", &EQ_OneParameter::ReadModelParameter));
  fKeyword.insert(pair<string, EQ_OP_MthPtr>(
      "k_targetparameterstdev", &EQ_OneParameter::ReadTargetParameterStDev));
  fKeyword.insert(pair<string, EQ_OP_MthPtr>(
      "k_nonzaifortmva", &EQ_OneParameter::ReadNonZaiTMVAVariables));
  DBGL;
}
//________________________________________________________________________
void EQ_OneParameter::ReadTargetParameter(const string &line) {
  DBGL;
  int pos = 0;
  string keyword = tlc(StringLine::NextWord(line, pos, ' '));
  if (keyword != "k_targetparameter")  // Check the keyword
  {
    ERROR << " Bad keyword : \"k_targetparameter\" not found !" << endl;
    exit(1);
  }

  fTargetParameter = StringLine::NextWord(line, pos, ' ');

  DBGL;
}

//________________________________________________________________________
void EQ_OneParameter::ReadModelParameter(const string &line) {
  DBGL;

  int pos = 0;
  string keyword = tlc(StringLine::NextWord(line, pos, ' '));
  if (keyword != "k_modelparameter")  // Check the keyword
  {
    ERROR << " Bad keyword : \"k_modelparameter\" not found !" << endl;
    exit(1);
  }

  keyword = StringLine::NextWord(line, pos, ' ');

  fModelParameter[keyword] = -1;

  DBGL;
}
//________________________________________________________________________
void EQ_OneParameter::ReadNonZaiTMVAVariables(const string &line) {
  DBGL;

  int pos = 0;
  string keyword = tlc(StringLine::NextWord(line, pos, ' '));
  if (keyword != "k_nonzaifortmva")  // Check the keyword
  {
    ERROR << " Bad keyword : \"k_nonzaifortmva\" not found !" << endl;
    exit(1);
  }

  keyword = StringLine::NextWord(line, pos, ' ');
  fListOfNonZaiTMVAVariables.push_back(make_pair(-1.0, keyword));
  DBGL;
}
//________________________________________________________________________
void EQ_OneParameter::SetNonZaiTMVAVariable(string snZP, double dnZP) {
  DBGL;
  for (unsigned int j = 0; j < fListOfNonZaiTMVAVariables.size(); j++) {
    if (fListOfNonZaiTMVAVariables[j].second == snZP) {
      fListOfNonZaiTMVAVariables[j].first = dnZP;
      return;
    }
  }
  fListOfNonZaiTMVAVariables.push_back(make_pair(dnZP, snZP));

  DBGL;
}
//________________________________________________________________________
void EQ_OneParameter::ReadTargetParameterStDev(const string &line) {
  DBGL;

  int pos = 0;
  string keyword = tlc(StringLine::NextWord(line, pos, ' '));
  if (keyword != "k_targetparameterstdev")  // Check the keyword
  {
    ERROR << " Bad keyword : \"k_targetparameterstdev\" not found !" << endl;
    exit(1);
  }

  fTargetParameterStDev = atof(StringLine::NextWord(line, pos, ' ').c_str());
  ;

  DBGL;
}

//________________________________________________________________________
void EQ_OneParameter::PrintInfo() {
  INFO << "Reactor Type : " << fDBRType << endl;
  INFO << "Fuel Type : " << fDBFType << endl;
  INFO << "Specific Power [W/g]: " << fSpecificPower << endl;

  map<string, IsotopicVector>::iterator it_s_IV;
  map<string, double>::iterator it_s_D;

  for (it_s_IV = fStreamList.begin(); it_s_IV != fStreamList.end(); it_s_IV++) {
    INFO << (*it_s_IV).first << "  (Z A I) :" << endl;
    map<ZAI, double>::iterator it1;
    map<ZAI, double> fMap1 =
        fStreamList[(*it_s_IV).first].GetIsotopicQuantity();
    for (it1 = fMap1.begin(); it1 != fMap1.end(); it1++)
      INFO << (*it1).first.Z() << " " << (*it1).first.A() << " "
           << (*it1).first.I() << endl;
  }
  INFO << "Minimum fraction in the fuel for each material : " << endl;
  for (it_s_D = fStreamListEqMMassFractionMin.begin();
       it_s_D != fStreamListEqMMassFractionMin.end(); it_s_D++) {
    INFO << (*it_s_D).first << " "
         << fStreamListEqMMassFractionMin[(*it_s_D).first] << endl;
  }

  INFO << "Maximum fraction in the fuel for each material : " << endl;
  for (it_s_D = fStreamListEqMMassFractionMax.begin();
       it_s_D != fStreamListEqMMassFractionMax.end(); it_s_D++) {
    INFO << (*it_s_D).first << " "
         << fStreamListEqMMassFractionMax[(*it_s_D).first] << endl;
  }

  INFO << "ZAI limits (validity domain)[prop in fresh fuel] (Z A I min max) :"
       << endl;
  for (map<ZAI, pair<double, double> >::iterator Domain_it = fZAILimits.begin();
       Domain_it != fZAILimits.end(); Domain_it++) {
    double ThatZAIMin = Domain_it->second.first;
    double ThatZAIMax = Domain_it->second.second;
    int Z = Domain_it->first.Z();
    int A = Domain_it->first.A();
    int I = Domain_it->first.I();

    INFO << ThatZAIMin << " < ZAI (" << Z << " " << A << " " << I << ")"
         << " < " << ThatZAIMax << endl;
  }
}
