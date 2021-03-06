#include "SeparationPlant.hxx"

#include "CLASSLogger.hxx"
#include "Scenario.hxx"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

//________________________________________________________________________
//
//		SeparationPlant
//
//
//
//
//________________________________________________________________________
ClassImp(SeparationPlant)

    SeparationPlant::SeparationPlant()
    : CLASSBackEnd(-2) {
  fOutBackEndFacility = 0;
  SetName("C_SeparationPlant.");
  SetIsStorageType();
}

//________________________________________________________________________
SeparationPlant::SeparationPlant(CLASSLogger* log) : CLASSBackEnd(log, -2) {
  fCycleTime = 0;
  fIsStarted = false;
  fPutToWaste = true;

  fOutBackEndFacility = 0;
  SetName("C_SeparationPlant.");

  SetIsStorageType();

  INFO << " A new SeparationPlant has been define :" << endl;
  INFO << "\t The Separation Time set at\t " << (double)(fCycleTime / cYear)
       << " year" << endl;
  WARNING << " All Separated Fuel go directly to WASTE after cooling !! "
          << endl;
}

//________________________________________________________________________
SeparationPlant::~SeparationPlant() {}

//________________________________________________________________________

//________________________________________________________________________
void SeparationPlant::SetBackEndDestination(CLASSBackEnd* storagedestination,
                                            IsotopicVector isotopicvector,
                                            cSecond destinationstartingtime) {
  DBGL;

  fDestinationStorageStartingTime.push_back(destinationstartingtime);
  fDestinationStorage.push_back(storagedestination);
  fDestinationStorageIV.push_back(isotopicvector);

  if (fDestinationStorage.size() != fDestinationStorageIV.size())
    ERROR << " fDestinationStorage.size() != fDestinationStorageIV.size() !! "
          << endl;

  DBGL;
}

//________________________________________________________________________
void SeparationPlant::AddIV(IsotopicVector IV) {
  DBGL;
  for (int fds = 0; fds < (int)fDestinationStorage.size(); fds++) {
    cSecond CurrentTime = GetParc()->GetAbsoluteTime();

    DBGV("Separation..." << endl);
    DBGV("Current Time : " << CurrentTime << endl);
    DBGV("IV Separation Time : " << fDestinationStorageStartingTime[fds]
                                 << endl);

    if (CurrentTime >= fDestinationStorageStartingTime[fds]) {
      IsotopicVector IVtmp;
      IVtmp = IV * fDestinationStorageIV[fds];
      fDestinationStorage[fds]->AddIV(IVtmp);
      IV -= IVtmp;
    }
  }

  GetParc()->AddWaste(IV);

  DBGL;
}

//________________________________________________________________________
map<cSecond, int> SeparationPlant::GetTheBackEndTimePath() {
  DBGL;

  map<cSecond, int> TheBackEndTimePath;
  for (int i = 0; i < (int)fDestinationStorage.size(); i++) {
    map<cSecond, int> TheBackEndTimePath_tmp =
        fDestinationStorage[i]->GetTheBackEndTimePath();
    map<cSecond, int>::iterator it;
    for (it = TheBackEndTimePath_tmp.begin();
         it != TheBackEndTimePath_tmp.end(); it++) {
      pair<map<cSecond, int>::iterator, bool> IResult;
      IResult = TheBackEndTimePath.insert(
          pair<cSecond, int>((*it).first, (*it).second));
      if (!IResult.second) IResult.first->second |= (*it).second;
    }
  }

  DBGL;
  return TheBackEndTimePath;
}

//________________________________________________________________________
