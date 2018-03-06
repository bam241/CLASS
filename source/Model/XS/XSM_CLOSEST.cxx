#include "XSM_CLOSEST.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"
#include "XSModel.hxx"

#include <TGraph.h>
#include <TString.h>
#include "TROOT.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

//________________________________________________________________________
//
//		XSM_CLOSEST
//
//
//
//
//________________________________________________________________________

XSM_CLOSEST::XSM_CLOSEST(string DB_index_file, bool oldreadmethod)
    : XSModel(new CLASSLogger("XSM_CLOSEST.log")) {
  fOldReadMethod = oldreadmethod;
  fDataBaseIndex = DB_index_file;
  fDistanceType = 0;
  fWeightedDistance = false;
  fEvolutionDataInterpolation = false;
  ReadDataBase();

  // Warning
  INFO << " A EvolutionData has been define :" << endl;
  INFO << "\t His index is : \"" << DB_index_file << "\" " << endl;
  INFO << "\t " << fFuelDataBank.size() << " EvolutionData have been read."
       << endl
       << endl;
}

//________________________________________________________________________
XSM_CLOSEST::XSM_CLOSEST(CLASSLogger* Log, string DB_index_file,
                         bool oldreadmethod)
    : XSModel(Log) {
  fOldReadMethod = oldreadmethod;
  fDataBaseIndex = DB_index_file;
  fDistanceType = 0;
  fWeightedDistance = false;
  fEvolutionDataInterpolation = false;
  ReadDataBase();

  // Warning
  INFO << " A EvolutionData has been define :" << endl;
  INFO << "\t His index is : \"" << DB_index_file << "\" " << endl;
  INFO << "\t " << fFuelDataBank.size() << " EvolutionData have been read."
       << endl
       << endl;
}

//________________________________________________________________________
XSM_CLOSEST::~XSM_CLOSEST() {
  for (int i = 0; i < (int)fFuelDataBank.size(); i++)
    fFuelDataBank[i].DeleteEvolutionData();
  fFuelDataBank.clear();
}

//________________________________________________________________________
void XSM_CLOSEST::ReadDataBase() {
  if (fFuelDataBank.size() != 0) {
    for (int i = 0; i < (int)fFuelDataBank.size(); i++)
      fFuelDataBank[i].DeleteEvolutionData();
    fFuelDataBank.clear();
  }

  ifstream DataDB(fDataBaseIndex.c_str());  // Open the File
  if (!DataDB) WARNING << " Can't open \"" << fDataBaseIndex << "\"\n" << endl;

  vector<double> vTime;
  vector<double> vTimeErr;

  string line;
  int start = 0;

  // First Get Fuel Type
  getline(DataDB, line);
  if (StringLine::NextWord(line, start, ' ') != "TYPE") {
    ERROR << " Bad Database file : " << fDataBaseIndex
          << " Can't find the type of the DataBase" << endl;
    exit(1);
  }
  fFuelType = StringLine::NextWord(line, start, ' ');

  // Then Get All the Database

  while (!DataDB.eof()) {
    getline(DataDB, line);
    if (line != "") {
      EvolutionData evolutionproduct(GetLog(), line, fOldReadMethod);
      fFuelDataBank.push_back(evolutionproduct);
    }
  }
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
/*			Distance Calculation			 */
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
map<double, int> XSM_CLOSEST::GetDistancesTo(IsotopicVector isotopicvector,
                                             double t) {
  map<double, int> distances;

  for (int i = 0; i < (int)fFuelDataBank.size(); i++) {
    pair<map<double, int>::iterator, bool> IResult;

    IsotopicVector ActinidesCompositionAtT =
        fFuelDataBank[i].GetIsotopicVectorAt(t).GetActinidesComposition();
    IsotopicVector IV_ActinidesComposition =
        isotopicvector.GetActinidesComposition();

    double NormalisationFactor =
        Norme(IV_ActinidesComposition) / Norme(ActinidesCompositionAtT);

    double distance = Distance(IV_ActinidesComposition,
                               ActinidesCompositionAtT / NormalisationFactor,
                               fDistanceType, fDistanceParameter);

    IResult = distances.insert(pair<double, int>(distance, i));
  }

  return distances;
}
//________________________________________________________________________
EvolutionData XSM_CLOSEST::GetCrossSections(IsotopicVector isotopicvector,
                                            double t) {
  DBGL;
  double distance = 0;
  int N_BestEvolutionData = 0;

  if (fWeightedDistance) {
    DBGL;
    IsotopicVector ActinidesCompositionAtT =
        fFuelDataBank[0].GetIsotopicVectorAt(t).GetActinidesComposition();
    IsotopicVector IV_ActinidesComposition =
        isotopicvector.GetActinidesComposition();

    double NormalisationFactor =
        Norme(ActinidesCompositionAtT) / Norme(IV_ActinidesComposition);

    distance = Distance(IV_ActinidesComposition / NormalisationFactor,
                        fFuelDataBank[0]);

    for (int i = 1; i < (int)fFuelDataBank.size(); i++) {
      double D = 0;
      ActinidesCompositionAtT =
          fFuelDataBank[i].GetIsotopicVectorAt(t).GetActinidesComposition();
      IV_ActinidesComposition = isotopicvector.GetActinidesComposition();
      NormalisationFactor =
          Norme(ActinidesCompositionAtT) / Norme(IV_ActinidesComposition);

      D = Distance(IV_ActinidesComposition / NormalisationFactor,
                   fFuelDataBank[i]);

      if (D < distance) {
        distance = D;
        N_BestEvolutionData = i;
      }
    }
    DBGL;
    return fFuelDataBank[N_BestEvolutionData];
  } else if (fEvolutionDataInterpolation) {
    DBGL;
    map<double, int> distance_map = GetDistancesTo(isotopicvector, t);
    map<double, int>::iterator it_distance;
    int NClose = 64;
    int Nstep = 0;
    EvolutionData EvolInterpolate;
    double SumOfDistance = 0;
    for (it_distance = distance_map.begin(); it_distance != distance_map.end();
         it_distance++) {
      double distance = (*it_distance).first;
      int ED_Indice = (*it_distance).second;

      if (distance == 0) {
        EvolInterpolate = Multiply(1, fFuelDataBank[ED_Indice]);
        return EvolInterpolate;
      }

      if (Nstep == 0)
        EvolInterpolate = Multiply(1. / distance, fFuelDataBank[ED_Indice]);
      else {
        EvolutionData Evoltmp = EvolInterpolate;
        EvolutionData Evoltmp2 =
            Multiply(1. / distance, fFuelDataBank[ED_Indice]);

        EvolInterpolate = Sum(Evoltmp, Evoltmp2);
        Evoltmp.DeleteEvolutionData();
        Evoltmp2.DeleteEvolutionData();
      }

      SumOfDistance += 1. / distance;
      Nstep++;
      if (Nstep == NClose) break;
    }

    EvolutionData Evoltmp = EvolInterpolate;
    EvolInterpolate.Clear();

    EvolInterpolate = Multiply(1 / SumOfDistance, Evoltmp);

    Evoltmp.DeleteEvolutionData();
    DBGL;
    return EvolInterpolate;

  } else {
    DBGL;
    IsotopicVector ActinidesCompositionAtT =
        fFuelDataBank[0].GetIsotopicVectorAt(t).GetActinidesComposition();
    IsotopicVector IV_ActinidesComposition =
        isotopicvector.GetActinidesComposition();

    double NormalisationFactor =
        Norme(IV_ActinidesComposition) / Norme(ActinidesCompositionAtT);

    distance = Distance(IV_ActinidesComposition,
                        ActinidesCompositionAtT / NormalisationFactor,
                        fDistanceType, fDistanceParameter);

    for (int i = 1; i < (int)fFuelDataBank.size(); i++) {
      double D = 0;
      ActinidesCompositionAtT =
          fFuelDataBank[i].GetIsotopicVectorAt(t).GetActinidesComposition();
      IV_ActinidesComposition = isotopicvector.GetActinidesComposition();

      NormalisationFactor =
          Norme(IV_ActinidesComposition) / Norme(ActinidesCompositionAtT);
      D = Distance(IV_ActinidesComposition,
                   ActinidesCompositionAtT / NormalisationFactor, fDistanceType,
                   fDistanceParameter);

      if (D < distance) {
        distance = D;
        N_BestEvolutionData = i;
      }
    }

    INFO << distance << endl;
    DBGL;

    return fFuelDataBank[N_BestEvolutionData];
  }
}
//________________________________________________________________________
void XSM_CLOSEST::CalculateDistanceParameter() {
  DBGL;
  if (fDistanceType != 1) {
    WARNING << " Distance Parameter will be calculate even if the distance "
               "type is not the good one. Any Distance Parameters given by the "
               "user will be overwriten"
            << endl;
  }

  fDistanceParameter.Clear();

  // We calculate the weight for the distance calculation.
  int NevolutionDatainFuelDataBank = 0;

  for (int i = 0; i < (int)fFuelDataBank.size(); i++) {
    NevolutionDatainFuelDataBank++;
    map<ZAI, double>::iterator itit;
    map<ZAI, double> isovector =
        fFuelDataBank[i].GetIsotopicVectorAt(0).GetIsotopicQuantity();
    for (itit = isovector.begin(); itit != isovector.end();
         itit++)  // Boucle sur ZAI
    {
      double TmpXS = 0;

      for (int i = 1; i < 4;
           i++)  // Loop on Reactions 1 == fission, 2 == capture, 3 == n2n
        TmpXS += fFuelDataBank[i].GetXSForAt(0, (*itit).first, i);

      fDistanceParameter.Add((*itit).first, TmpXS);
    }
  }
  fDistanceParameter.Multiply((double)1.0 / NevolutionDatainFuelDataBank);

  if (GetLog()) {
    INFO << "!!INFO!! Distance Parameters " << endl;
    map<ZAI, double>::iterator it2;
    for (it2 = fDistanceParameter.GetIsotopicQuantity().begin();
         it2 != fDistanceParameter.GetIsotopicQuantity().end(); it2++) {
      INFO << (*it2).first.Z() << " ";
      INFO << (*it2).first.A() << " ";
      INFO << (*it2).first.I() << " ";
      INFO << ": " << (*it2).second;
      INFO << endl;
    }
    INFO << endl;
  }

  DBGL;
}
//________________________________________________________________________
void XSM_CLOSEST::SetDistanceParameter(IsotopicVector DistanceParameter) {
  fDistanceParameter = DistanceParameter;

  INFO << "!!INFO!! Distance Parameters " << endl;
  map<ZAI, double>::iterator it2;
  for (it2 = fDistanceParameter.GetIsotopicQuantity().begin();
       it2 != fDistanceParameter.GetIsotopicQuantity().end(); it2++) {
    INFO << (*it2).first.Z() << " ";
    INFO << (*it2).first.A() << " ";
    INFO << (*it2).first.I() << " ";
    INFO << ": " << (*it2).second;
    INFO << endl;
  }
  INFO << endl;
}

//________________________________________________________________________
void XSM_CLOSEST::SetDistanceType(int DistanceType) {
  fDistanceType = DistanceType;
  if (fDistanceType == 1) {
    CalculateDistanceParameter();
  } else if (fDistanceType == 2 && Norme(fDistanceParameter) == 0) {
    // This is so bad!! You will probably unsynchronize all the reactor....
    WARNING << " Distance use weight defined by user for each isotope, but no "
               "weight have been given"
            << endl
            << "Use SetDistanceParameter()" << endl;
    exit(1);
  } else if (fDistanceType != 0 && fDistanceType != 1 && fDistanceType != 2) {
    ERROR << " Distancetype defined by the user isn't recognized by the code"
          << endl;
    exit(1);
  }
}
