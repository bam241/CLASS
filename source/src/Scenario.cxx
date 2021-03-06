#include "Scenario.hxx"
#include "CLASSMethod.hxx"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "stdlib.h"
#include "time.h"

#include "CLASSBackEnd.hxx"
#include "CLASSLogger.hxx"
#include "FabricationPlant.hxx"
#include "Pool.hxx"
#include "Reactor.hxx"
#include "SeparationPlant.hxx"
#include "Storage.hxx"

//________________________________________________________________________
//
//		CLASS
//
//
//
//
//________________________________________________________________________

//________________________________________________________________________
Scenario::Scenario(CLASSLogger* log, cSecond abstime) : CLASSObject(log) {
  fNewTtree = true;
  fPrintStep = (cSecond)(cYear);  // One Step per Year
  fAbsoluteTime = abstime;
  fStartingTime = fAbsoluteTime;

  fStockManagement = true;
  fLogTimeStep = false;

  fOutputFileName = "CLASS_Default.root";
  fOutputTreeName = "Data";
  fOutFile = 0;
  fOutT = 0;

  fParcPower = 0;

  fZAIThreshold = -1;

  // Warning

  INFO << " Parc has been define :" << endl;
  INFO << " Print  set at : " << (double)(fPrintStep / cYear) << " year"
       << endl;
  INFO << " StockManagement set at : true" << endl;
  INFO << " OutPut will be in \"" << fOutputFileName << "\" File and \""
       << fOutputTreeName << "\" TTree" << endl;
  INFO << " Log will be in " << GetLog()->GetCLASSLoggerName() << endl;
}
//________________________________________________________________________
Scenario::Scenario(cSecond abstime, CLASSLogger* log) : CLASSObject(log) {
  fNewTtree = true;
  fPrintStep = (cSecond)(cYear);  // One Step per Year
  fAbsoluteTime = abstime;
  fStartingTime = fAbsoluteTime;

  fStockManagement = true;
  fLogTimeStep = false;

  fOutputFileName = "CLASS_Default.root";
  fOutputTreeName = "Data";
  fOutFile = 0;
  fOutT = 0;

  fParcPower = 0;

  fZAIThreshold = -1;

  // Warning

  INFO << " Parc has been define :" << endl;
  INFO << " Print  set at : " << (double)(fPrintStep / cYear) << " year"
       << endl;
  INFO << " StockManagement set at : true" << endl;
  INFO << " OutPut will be in \"" << fOutputFileName << "\" File and \""
       << fOutputTreeName << "\" TTree" << endl;
  INFO << " Log will be in " << GetLog()->GetCLASSLoggerName() << endl;
}
//________________________________________________________________________
Scenario::Scenario(cSecond abstime)
    : CLASSObject(new CLASSLogger("CLASS_OUTPUT.log")) {
  fNewTtree = true;
  fPrintStep = (cSecond)(cYear);  // One Step per Year
  fAbsoluteTime = abstime;
  fStartingTime = fAbsoluteTime;

  fStockManagement = true;
  fLogTimeStep = false;

  fOutputFileName = "CLASS_Default.root";
  fOutputTreeName = "Data";
  fOutFile = 0;
  fOutT = 0;

  fParcPower = 0;

  fZAIThreshold = -1;

  // Warning

  INFO << " Parc has been define :" << endl;
  INFO << " Print  set at : " << (double)(fPrintStep / cYear) << " year"
       << endl;
  INFO << " StockManagement set at : true" << endl;
  INFO << " OutPut will be in \"" << fOutputFileName << "\" File and \""
       << fOutputTreeName << "\" TTree" << endl;
  INFO << " Log will be in " << GetLog()->GetCLASSLoggerName() << endl;
}
//________________________________________________________________________
Scenario::~Scenario() {
#pragma omp single
  { CloseOutputTree(); }
}
//________________________________________________________________________
void Scenario::AddPool(Pool* Pool) {
  fPool.push_back(Pool);
  fPool.back()->SetParc(this);
  fPool.back()->SetDecayDataBank((*this).GetDecayDataBase());
  fPool.back()->SetId((int)fPool.size() - 1);
  fPool.back()->SetInternalTime(fAbsoluteTime);

  string Pool_name = fPool.back()->GetName();
  if (Pool_name == "P_Pool.") {
    Pool_name = "P_Pool";
    Pool_name += dtoa(fPool.back()->GetId());
    Pool_name += ".";
    fPool.back()->SetName(Pool_name.c_str());
  } else {
    string name_tmp = Pool_name;
    Pool_name = "P_";
    Pool_name += name_tmp;
    Pool_name += ".";
    fPool.back()->SetName(Pool_name.c_str());
  }

  if (!fNewTtree) fOutT->Branch(fPool.back()->GetName(), "Pool", &fPool.back());
}
//________________________________________________________________________
void Scenario::AddReactor(Reactor* reactor) {
  fReactor.push_back(reactor);
  fReactor.back()->SetParc(this);
  fReactor.back()->SetId((int)fReactor.size() - 1);
  fReactor.back()->SetInternalTime(fAbsoluteTime);

  string Reactor_name = fReactor.back()->GetName();
  if (Reactor_name == "R_Reactor.") {
    Reactor_name = "R_Reactor";
    Reactor_name += dtoa(fReactor.back()->GetId());
    Reactor_name += ".";
    fReactor.back()->SetName(Reactor_name.c_str());
  } else {
    string name_tmp = Reactor_name;
    Reactor_name = "R_";
    Reactor_name += name_tmp;
    Reactor_name += ".";
    fReactor.back()->SetName(Reactor_name.c_str());
  }

  if (!fNewTtree)
    fOutT->Branch(fReactor.back()->GetName(), "Reactor", &fReactor.back());
}
//________________________________________________________________________
void Scenario::AddStorage(Storage* storage) {
  fStorage.push_back(storage);
  fStorage.back()->SetParc(this);
  fStorage.back()->SetDecayDataBank((*this).GetDecayDataBase());
  fStorage.back()->SetId((int)fStorage.size() - 1);
  fStorage.back()->SetInternalTime(fAbsoluteTime);

  string Storage_name = fStorage.back()->GetName();

  if (Storage_name == "S_Storage.") {
    Storage_name = "S_Storage";
    Storage_name += dtoa(fStorage.back()->GetId());
    Storage_name += ".";
    fStorage.back()->SetName(Storage_name.c_str());
  } else {
    string name_tmp = Storage_name;
    Storage_name = "S_";
    Storage_name += name_tmp;
    Storage_name += ".";
    fStorage.back()->SetName(Storage_name.c_str());
  }

  if (!fNewTtree)
    fOutT->Branch(fStorage.back()->GetName(), "Storage", &fStorage.back());
}
//________________________________________________________________________
void Scenario::AddFabricationPlant(FabricationPlant* fabricationplant) {
  fFabricationPlant.push_back(fabricationplant);
  fFabricationPlant.back()->SetParc(this);
  fFabricationPlant.back()->SetDecayDataBank((*this).GetDecayDataBase());
  fFabricationPlant.back()->SetId((int)fStorage.size() - 1);
  fFabricationPlant.back()->SetInternalTime(fAbsoluteTime);

  string FP_name = fFabricationPlant.back()->GetName();
  if (FP_name == "F_FabricationPlant.") {
    FP_name = "F_FabricationPlant";
    FP_name += dtoa(fFabricationPlant.back()->GetId());
    FP_name += ".";
    fFabricationPlant.back()->SetName(FP_name.c_str());
  } else {
    string name_tmp = FP_name;
    FP_name = "F_";
    FP_name += name_tmp;
    FP_name += ".";
    fFabricationPlant.back()->SetName(FP_name.c_str());
  }

  if (!fNewTtree)
    fOutT->Branch(fFabricationPlant.back()->GetName(), "FabricationPlant",
                  &fFabricationPlant.back());
}
//________________________________________________________________________
void Scenario::AddSeparationPlant(SeparationPlant* SeparationPlant) {
  fSeparationPlant.push_back(SeparationPlant);
  fSeparationPlant.back()->SetParc(this);
  fSeparationPlant.back()->SetDecayDataBank((*this).GetDecayDataBase());
  fSeparationPlant.back()->SetId((int)fSeparationPlant.size() - 1);
  fSeparationPlant.back()->SetInternalTime(fAbsoluteTime);

  string SeparationPlant_name = fSeparationPlant.back()->GetName();
  if (SeparationPlant_name == "C_SepPlant.") {
    SeparationPlant_name = "C_SepPlant";
    SeparationPlant_name += dtoa(fSeparationPlant.back()->GetId());
    SeparationPlant_name += ".";
    fSeparationPlant.back()->SetName(SeparationPlant_name.c_str());
  } else {
    string name_tmp = SeparationPlant_name;
    SeparationPlant_name = "C_";
    SeparationPlant_name += name_tmp;
    SeparationPlant_name += ".";
    fSeparationPlant.back()->SetName(SeparationPlant_name.c_str());
  }

  if (!fNewTtree)
    fOutT->Branch(fSeparationPlant.back()->GetName(), "SeparationPlant",
                  &fSeparationPlant.back());
}

//________________________________________________________________________

//________________________________________________________________________
void Scenario::BuildTimeVector(cSecond t) {
  DBGL;

  fTimeStep.clear();
  fTimeStep.insert(pair<cSecond, int>(t, 1));
  //********* Printing Step *********//
  {
    DBGL;
    cSecond step = fStartingTime;

    if (step >= fAbsoluteTime) fTimeStep.insert(pair<cSecond, int>(step, 1));
    step += fPrintStep;
    cSecond timescale = 1;

    do {
      if (step >= fAbsoluteTime) fTimeStep.insert(pair<cSecond, int>(step, 1));

      if (fLogTimeStep) {
        timescale *= 10;
        step = fPrintStep * timescale;
      } else
        step += fPrintStep;

    } while (step < t);
    DBGL;
  }

  for (int i = 0; i < (int)fReactor.size(); i++) {
    DBGL;
    cSecond R_StartingTime = fReactor[i]->GetCreationTime();
    cSecond R_ShutDownTime =
        fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime();

    ScheduleEntry* R_Entry =
        fReactor[i]->GetScheduler()->GetEntryAt(R_StartingTime);
    double R_BU = R_Entry->GetBurnUp();
    double R_Power = R_Entry->GetPower();
    double R_HMMass = R_Entry->GetHeavyMetalMass();

    cSecond R_CycleTime =
        (cSecond)(R_BU * 1e9 / R_Power * R_HMMass * 3600 * 24);

    DBGV("BU : " << R_BU << " R_Power : " << R_Power << " R_R_HMMass : "
                 << R_HMMass << " R_CycleTime : " << R_CycleTime);

    if (R_CycleTime == 0) {
      ERROR << " Be carefull a reactor cycletime is set to 0 second....\"\n"
            << endl;
      exit(1);
    }

    int R_FacilityType = fReactor[i]->GetFacilityType();

    cSecond F_CycleTime = 0;

    cSecond step = R_StartingTime;

    map<cSecond, int> R_BackEndTimePath =
        fReactor[i]->GetOutBackEndFacility()->GetTheBackEndTimePath();

    if (R_Entry->GetReactorModel()->GetPhysicsModels())
      F_CycleTime = fReactor[i]->GetFabricationPlant()->GetCycleTime();

    //********* Reactor Evolution Step *********//
    // ShutDown of a reactor

    // Test if the shutdown of the reactor is after the actual time
    // (AbsolutreTime) and before the end of the evolution (t)
    if (R_ShutDownTime < t) {
      // Shutdown
      if (R_ShutDownTime > fAbsoluteTime) {
        pair<map<cSecond, int>::iterator, bool> IResult;
        IResult = fTimeStep.insert(pair<cSecond, int>(R_ShutDownTime, 2));
        if (!IResult.second) IResult.first->second |= 2;
      }

      // BackEnd fuel Cycle after reactor Shutdown

      map<cSecond, int>::iterator TV_it;  // the time vector iterator

      // Loop on the BackEnd fuel Cycle Time path
      for (TV_it = R_BackEndTimePath.begin(); TV_it != R_BackEndTimePath.end();
           TV_it++) {
        // Test if each step of the Fuel Cycle BackEnd is after the actual time
        // (AbsolutreTime) and before the end of the evolution (t)
        if (R_ShutDownTime + (*TV_it).first >= fAbsoluteTime &&
            R_ShutDownTime + (*TV_it).first <= t) {
          pair<map<cSecond, int>::iterator, bool> IResult;
          IResult = fTimeStep.insert(pair<cSecond, int>(
              R_ShutDownTime + (*TV_it).first, (*TV_it).second));
          if (!IResult.second) IResult.first->second |= (*TV_it).second;
        }
      }
    }

    // Start the reactor and the Fuel Fabrication
    if (step >= fAbsoluteTime && step <= t && step < R_ShutDownTime) {
      pair<map<cSecond, int>::iterator, bool> IResult;
      IResult = fTimeStep.insert(pair<cSecond, int>(step, R_FacilityType));
      if (!IResult.second) IResult.first->second |= R_FacilityType;
    }

    //********* FabricationPlant Evolution Step *********//

    if (R_Entry->GetReactorModel()->GetPhysicsModels()) {
      fReactor[i]->GetFabricationPlant()->AddReactor(i, step);

      F_CycleTime = fReactor[i]->GetFabricationPlant()->GetCycleTime();

      if (step - F_CycleTime >= fAbsoluteTime && step - F_CycleTime <= t &&
          step < R_ShutDownTime) {  // Set End of reactor cycle
        pair<map<cSecond, int>::iterator, bool> IResult;
        IResult = fTimeStep.insert(pair<cSecond, int>(step - F_CycleTime, 16));
        if (!IResult.second) IResult.first->second |= 16;
      } else if (step - F_CycleTime < fStartingTime) {
        ERROR << " Can't Build Fuel before Scenario's start\"\n" << endl;
        exit(1);
      }
    }

    step += R_CycleTime;
    // Prepare the first Cycle
    R_Entry = fReactor[i]->GetScheduler()->GetEntryAt(step);
    R_BU = R_Entry->GetBurnUp();
    R_Power = R_Entry->GetPower();
    R_HMMass = R_Entry->GetHeavyMetalMass();
    R_CycleTime = (cSecond)(R_BU * 1e9 / R_Power * R_HMMass * 3600 * 24);

    DBGV("BU : " << R_BU << " R_Power : " << R_Power << " R_R_HMMass : "
                 << R_HMMass << " R_CycleTime : " << R_CycleTime);
    if (R_CycleTime == 0) {
      ERROR << " Be carefull a reactor cycletime is set to 0 second....\"\n"
            << endl;
      exit(1);
    }

    while (step <= t && step <= R_ShutDownTime) {
      DBGL;

      // FabricationPlant Evolution Step
      if (R_Entry->GetReactorModel()->GetPhysicsModels()) {
        fReactor[i]->GetFabricationPlant()->AddReactor(i, step);

        F_CycleTime = fReactor[i]->GetFabricationPlant()->GetCycleTime();

        if (step - F_CycleTime >= fAbsoluteTime && step - F_CycleTime <= t &&
            step < R_ShutDownTime) {  // Set End of reactor cycle
          pair<map<cSecond, int>::iterator, bool> IResult;
          IResult =
              fTimeStep.insert(pair<cSecond, int>(step - F_CycleTime, 16));
          if (!IResult.second) IResult.first->second |= 16;
        }
      }

      if (step >= fAbsoluteTime && step <= t &&
          step < R_ShutDownTime) {  // Set End of reactor cycle
        pair<map<cSecond, int>::iterator, bool> IResult =
            fTimeStep.insert(pair<cSecond, int>(step, 4));
        if (!IResult.second) IResult.first->second |= 4;
      }

      // End/Start Of Reactor Cycle Step //

      map<cSecond, int>::iterator TV_it;  // the time vector iterator
      // BackEnd fuel Cycle
      // Loop on the BackEnd fuel Cycle Time path
      for (TV_it = R_BackEndTimePath.begin(); TV_it != R_BackEndTimePath.end();
           TV_it++) {
        if (step + (*TV_it).first >= fAbsoluteTime &&
            step + (*TV_it).first <= t) {  // Test if each step of the Fuel
                                           // Cycle BackEnd is after the actual
                                           // time (AbsolutreTime) and before
                                           // the end of the evolution (t)
          pair<map<cSecond, int>::iterator, bool> IResult;

          IResult = fTimeStep.insert(
              pair<cSecond, int>(step + (*TV_it).first, (*TV_it).second));
          if (!IResult.second) IResult.first->second |= (*TV_it).second;
        }
      }

      step += R_CycleTime;

      // Update to the next fuel
      R_Entry = fReactor[i]->GetScheduler()->GetEntryAt(step);
      double R_BU = R_Entry->GetBurnUp();
      double R_Power = R_Entry->GetPower();
      double R_HMMass = R_Entry->GetHeavyMetalMass();
      R_CycleTime = (cSecond)(R_BU * 1e9 / R_Power * R_HMMass * 3600 * 24);

      if (R_CycleTime == 0) {
        ERROR << " Be carefull a reactor cycletime is set to 0 second....\"\n"
              << endl;
        exit(1);
      }

      DBGL;
    }

    DBGL;
  }
  //****** Print the Time Index ******//
  ofstream TimeStepfile("CLASS_TimeStep", ios_base::app);  // Open the File

  if (!TimeStepfile) WARNING << " Can't open \" CLASS_TimeStep \"\n" << endl;

  map<cSecond, int>::iterator it;
  for (it = fTimeStep.begin(); it != fTimeStep.end(); it++)
    TimeStepfile << (*it).first << " " << (*it).second << endl;

  DBGL;
}

//________________________________________________________________________
//___________________________ Evolution Method ___________________________
//________________________________________________________________________

void Scenario::BackEndEvolution() {
  DBGL;
  StorageEvolution();
  PoolEvolution();

  PoolDump();
  DBGL;
}
//________________________________________________________________________
void Scenario::PoolEvolution() {
  DBGL;
#pragma omp parallel for
  for (int i = 0; i < (int)fPool.size(); i++)
    fPool[i]->Evolution(fAbsoluteTime);

  DBGL;
}
//________________________________________________________________________
void Scenario::StorageEvolution() {
  DBGL;
#pragma omp parallel for
  for (int i = 0; i < (int)fStorage.size(); i++)
    fStorage[i]->Evolution(fAbsoluteTime);

  DBGL;
}
//________________________________________________________________________
void Scenario::FabricationPlantEvolution() {
  DBGL;
  //#pragma omp parallel for
  for (int i = 0; i < (int)fFabricationPlant.size(); i++)
    fFabricationPlant[i]->Evolution(fAbsoluteTime);

  DBGL;
}
//________________________________________________________________________
void Scenario::PoolDump() {
  DBGL;
  for (int i = 0; i < (int)fPool.size(); i++) fPool[i]->Dump();
  DBGL;
}

//________________________________________________________________________
void Scenario::ReactorEvolution() {
  DBGL;
  fParcPower = 0;
  fParcElectricPower = 0;
#pragma omp parallel for
  for (int i = 0; i < (int)fReactor.size(); i++)
    fReactor[i]->Evolution(fAbsoluteTime);

  for (int i = 0; i < (int)fReactor.size(); i++) fReactor[i]->Dump();

  DBGL;
}

//________________________________________________________________________
void Scenario::Evolution(cSecond t) {
  DBGL;
  fCloverCount = 0;
  PrintCLASSPresentation();

  BuildTimeVector(t);

  if (fNewTtree) {
    OpenOutputTree();
    OutAttach();
    UpdateParc();
    fOutT->Fill();
  }

  map<cSecond, int>::iterator it;
  for (it = fTimeStep.begin(); it != fTimeStep.end(); it++) {
    fWaste = fDecayDataBase->GetDecay(fWaste, (*it).first - fAbsoluteTime);
    fAbsoluteTime = (*it).first;

    if ((*it).second & 1 || (*it).second & 2 || (*it).second & 4 ||
        (*it).second & 8 || (*it).second & 16)
      BackEndEvolution();

    if ((*it).second & 1 || (*it).second & 2 || (*it).second & 4 ||
        (*it).second & 16)
      FabricationPlantEvolution();

    if ((*it).second & 1 || (*it).second & 2 || (*it).second & 4)
      ReactorEvolution();

    if ((*it).second & 1 || it == fTimeStep.begin()) {
#pragma omp single
      {
        if (fZAIThreshold != -1) ApplyZAIThreshold();

        UpdateParc();
        fOutT->Fill();
        ProgressPrintout((cSecond)t);
      }
    }
  }
  cout << "Calculation complete successfully     " << endl << endl;

  DBGL;
}

//_______________________________________________________________________
//______________________________ Out Method ______________________________
//________________________________________________________________________

//________________________________________________________________________
void Scenario::ProgressPrintout(cSecond t) {
  double Time = (fAbsoluteTime - fStartingTime) / (double)cYear;
  double Total = (t - fStartingTime) / (double)cYear;

  // Reset the line
  for (int i = 0; i < 10; i++) cout << "  ";
  cout << flush;

  cout << "\r[";
  for (int i = 0; i < (int)(Time / Total * 20.0); i++) cout << "|";
  for (int i = 20; i >= (int)(Time / Total * 20.0); i--) cout << "-";
  cout << "] ";

  // cout << " Processed ";
  if (Time < 10) cout << " ";
  if (Time < 100) cout << " ";
  cout << (int)Time << " / " << (int)Total << " Years \r";
  if (fLog->GetVerboseLVL() < 2)
    cout << flush;
  else
    cout << endl;
}
//________________________________________________________________________
void Scenario::PrintCLASSPresentation() {
  cout << endl;
  cout << "╭────────────────────────────────────────────────╮" << endl;
  cout << "│   ██████╗██╗      █████╗ ███████╗███████╗       │" << endl;
  cout << "│  ██╔════╝██║     ██╔══██╗██╔════╝██╔════╝       │" << endl;
  cout << "│  ██║     ██║     ███████║███████╗███████╗       │" << endl;
  cout << "│  ██║     ██║     ██╔══██║╚════██║╚════██║       │" << endl;
  cout << "│  ╚██████╗███████╗██║  ██║███████║███████║       │" << endl;
  cout << "│   ╚═════╝╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝       │" << endl;
  cout << "│                                  Version 5.1    │" << endl;
  cout << "├───────────────────────────────────────────────  ┤" << endl;
  cout << "│ Core Library for Advanced Scenarios Simulations │" << endl;
  cout << "│                                                 │" << endl;
  cout << "│      https://gitlab.in2p3.fr/sens/CLASS         │" << endl;
  cout << "├─────────────────────────────────────────────────┤" << endl;
  cout << "│ CLASS is governed by the CeCILL-B license       │" << endl;
  cout << "│ See http://www.cecill.info                      │" << endl;
  cout << "╰────────────────────────────────────────────────╯" << endl;
  cout << "                                                   " << endl;
  cout << "        STARTING FUEL CYCLE EVOLUTION              " << endl;
  cout << "                                                   " << endl;
  cout << "Evolution progression   :                          " << endl;
}
//________________________________________________________________________
void Scenario::ApplyZAIThreshold() {
  for (int i = 0; i < (int)fFabricationPlant.size(); i++)
    fFabricationPlant[i]->ApplyZAIThreshold(fZAIThreshold);

  for (int i = 0; i < (int)fPool.size(); i++)
    fPool[i]->ApplyZAIThreshold(fZAIThreshold);

  for (int i = 0; i < (int)fStorage.size(); i++)
    fStorage[i]->ApplyZAIThreshold(fZAIThreshold);

  for (int i = 0; i < (int)fReactor.size(); i++)
    fReactor[i]->ApplyZAIThreshold(fZAIThreshold);

  fWaste.ApplyZAIThreshold(fZAIThreshold);
}
//________________________________________________________________________
void Scenario::UpdateParc() {
  ResetQuantity();

  for (int i = 0; i < (int)fFabricationPlant.size(); i++)
    fFuelFabrication += fFabricationPlant[i]->GetInsideIV();

  for (int i = 0; i < (int)fPool.size(); i++)
    fTotalCooling += fPool[i]->GetInsideIV();

  for (int i = 0; i < (int)fStorage.size(); i++)
    fTotalStorage += fStorage[i]->GetInsideIV();

  for (int i = 0; i < (int)fReactor.size(); i++)
    fTotalInReactor += fReactor[i]->GetIVReactor();

  fIVTotal = fWaste + fTotalStorage + fTotalCooling + fFuelFabrication +
             fTotalInReactor;
  fIVInCycleTotal =
      fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;
}
//________________________________________________________________________
void Scenario::ResetQuantity() {
  fTotalInReactor.Clear();
  fTotalStorage.Clear();
  fTotalCooling.Clear();
  fFuelFabrication.Clear();
  fIVInCycleTotal.Clear();
  fIVTotal.Clear();
}
//________________________________________________________________________
void Scenario::OpenOutputTree() {
  INFO << "Opening : " << fOutputFileName << " ...\t";
  fOutFile = new TFile(fOutputFileName.c_str(), "UPDATE");

  if (!fOutFile) {
    ERROR << "\nCould not open " << fOutputFileName << endl;
    exit(-1);
  }
  INFO << "\t ...OK!" << endl;

  fOutT = new TTree(fOutputTreeName.c_str(), "Data Tree");
  INFO << "Creating Data Tree ...\t";
  if (!fOutT) {
    ERROR << "\nCould not create Data Tree in " << fOutputFileName << endl;
    exit(-1);
  }
  fNewTtree = false;
  INFO << "\t ...OK!" << endl;
}
//________________________________________________________________________
void Scenario::CloseOutputTree() {
  fOutFile->ls();
  INFO << "Writing outTree " << fOutputFileName << endl;
  fOutFile->Write("", TObject::kOverwrite);

  if (fOutFile->IsOpen()) {
    INFO << "Deleting outTree : " << endl;
    delete fOutT;
    INFO << "Closing file : " << fOutputFileName << endl;
    fOutFile->Close();
    delete fOutFile;
  } else {
    ERROR << "File was not opened " << fOutputFileName << endl;
    exit(-1);
  }
}
//________________________________________________________________________
void Scenario::OutAttach() {
  ResetQuantity();
  // Branch Absolut Time
  fOutT->Branch("AbsTime", &fAbsoluteTime, "AbsoluteTime/L");
  // Branch The Power installed in the Parc
  fOutT->Branch("ParcPower", &fParcPower, "ParcPower/D");

  // Branch the Sum IV

  fOutT->Branch("STOCK.", "IsotopicVector", &fTotalStorage);
  fOutT->Branch("FUELFABRICATION.", "IsotopicVector", &fFuelFabrication);
  fOutT->Branch("COOLING.", "IsotopicVector", &fTotalCooling);
  fOutT->Branch("REACTOR.", "IsotopicVector", &fTotalInReactor);
  fOutT->Branch("INCYCLE.", "IsotopicVector", &fIVInCycleTotal);
  fOutT->Branch("TOTAL.", "IsotopicVector", &fIVTotal);

  fOutT->Branch("OUTINCOME.", "IsotopicVector", &fOutIncome);
  fOutT->Branch("WASTE.", "IsotopicVector", &fWaste);

  // Branch the separate object

  for (int i = 0; i < (int)fStorage.size(); i++)
    fOutT->Branch(fStorage[i]->GetName(), "Storage", &fStorage[i]);

  for (int i = 0; i < (int)fPool.size(); i++)

    fOutT->Branch(fPool[i]->GetName(), "Pool", &fPool[i]);

  for (int i = 0; i < (int)fReactor.size(); i++)

    fOutT->Branch(fReactor[i]->GetName(), "Reactor", &fReactor[i]);

  for (int i = 0; i < (int)fFabricationPlant.size(); i++)
    fOutT->Branch(fFabricationPlant[i]->GetName(), "FabricationPlant",
                  &fFabricationPlant[i]);
}
//________________________________________________________________________
void Scenario::Write() {}
//________________________________________________________________________
void Scenario::Print() {
  for (int i = 0; i < (int)fPool.size(); i++) {
    INFO << "!!!!!!!!!STEP : " << fAbsoluteTime / (int)(cYear) << endl;
    INFO << "Pool : " << endl;
    INFO << "Cooling ";
    INFO << fPool[i]->GetIVArray().size() << endl;
  }

  for (int i = 0; i < (int)fReactor.size(); i++) {
    INFO << "Reactor" << endl;
    fReactor[i]->GetIVReactor().Print();
  }
}
