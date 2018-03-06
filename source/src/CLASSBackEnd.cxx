#include "CLASSBackEnd.hxx"

#include "CLASSLogger.hxx"
#include "DecayDataBank.hxx"
#include "Scenario.hxx"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

//________________________________________________________________________
//
//		CLASSBackEnd
//
//
//
//
//________________________________________________________________________
ClassImp(CLASSBackEnd)

    CLASSBackEnd::CLASSBackEnd(int type)
    : CLASSFacility(type), fDecayDataBase(0) {
  ;
}

CLASSBackEnd::CLASSBackEnd(CLASSLogger* log, int type)
    : CLASSFacility(log, type), fDecayDataBase(0) {
  ;
}

CLASSBackEnd::CLASSBackEnd(CLASSLogger* log, cSecond cycletime, int type)
    : CLASSFacility(log, cycletime, type), fDecayDataBase(0) {
  ;
}

//________________________________________________________________________
void CLASSBackEnd::ClearIVArray() {
  fInsideIV = IsotopicVector();
  fIVArray.clear();
  fIVArrayArrivalTime.clear();
}

//________________________________________________________________________
void CLASSBackEnd::AddIV(IsotopicVector isotopicvector) {
  AddCumulativeIVIn(isotopicvector);

  fIVArray.push_back(isotopicvector);
  fIVArrayArrivalTime.push_back(fInternalTime);
}
//________________________________________________________________________
void CLASSBackEnd::UpdateInsideIV() {
  DBGL;
  fInsideIV = IsotopicVector();
  for (std::size_t i = 0; i < fIVArray.size(); ++i) {
    fInsideIV += fIVArray[i];
  }
  DBGL;
}

//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector CLASSBackEnd::GetDecay(IsotopicVector const& isotopicvector,
                                      cSecond t) {
  DBGL;

  IsotopicVector IV;

  for (IsotopicVector::const_iterator it = isotopicvector.begin();
       it != isotopicvector.end(); ++it) {
    if (it->second > 0) {
      IV += fDecayDataBase->Evolution(it->first, t) * (it->second);
    }
  }

  DBGL;
  return IV;
}

void CLASSBackEnd::ApplyZAIThreshold(int z) {
  fInsideIV.ApplyZAIThreshold(z);
  fCumulativeIVIn.ApplyZAIThreshold(z);
  fCumulativeIVOut.ApplyZAIThreshold(z);

  for (std::size_t i = 0; i < fIVArray.size(); ++i) {
    fIVArray[i].ApplyZAIThreshold(z);
  }
}

std::map<cSecond, int> CLASSBackEnd::GetTheBackEndTimePath() {
  DBGL;
  std::map<cSecond, int> TheBackEndTimePath;

  if (!fIsStorageType) {
    int FacilityType = GetFacilityType();
    cSecond step = GetCycleTime();

    std::pair<std::map<cSecond, int>::iterator, bool> IResult =
        TheBackEndTimePath.insert(std::pair<cSecond, int>(step, FacilityType));

    if (!IResult.second) {
      IResult.first->second |= FacilityType;
    }

    std::map<cSecond, int> TheBackEndTimePath_tmp =
        GetOutBackEndFacility()->GetTheBackEndTimePath();

    for (std::map<cSecond, int>::iterator it = TheBackEndTimePath_tmp.begin();
         it != TheBackEndTimePath_tmp.end(); ++it) {
      std::pair<std::map<cSecond, int>::iterator, bool> IResult;

      IResult = TheBackEndTimePath.insert(
          std::pair<cSecond, int>(step + it->first, it->second));

      if (!IResult.second) {
        IResult.first->second |= it->second;
      }
    }
  }

  DBGL;
  return TheBackEndTimePath;
}