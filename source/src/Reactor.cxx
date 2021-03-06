#include "Reactor.hxx"

#include "CLASSConstante.hxx"
#include "EvolutionData.hxx"
#include "FabricationPlant.hxx"
#include "Pool.hxx"
#include "Scenario.hxx"
#include "Storage.hxx"

#include <cmath>
#include <iostream>

#include <typeinfo>

//________________________________________________________________________
//
//		Reactor
//
//
//
//
//________________________________________________________________________

ClassImp(Reactor)

    Reactor::Reactor()
    : CLASSFacility(4) {
  SetName("R_Reactor.");

  fOutBackEndFacility = 0;
  fStorage = 0;
  fFabricationPlant = 0;
  fReactorScheduler = 0;
}
//________________________________________________________________________
Reactor::Reactor(CLASSLogger* log) : CLASSFacility(log, 4) {
  DBGL;

  fOutBackEndFacility = 0;
  fStorage = 0;
  fFabricationPlant = 0;
  fReactorScheduler = 0;
  SetName("R_Reactor.");

  DBGL;
}
//________________________________________________________________________
Reactor::Reactor(CLASSLogger* log, CLASSBackEnd* Pool, cSecond creationtime,
                 cSecond lifetime, double power, double HMMass,
                 double CapacityFactor)
    : CLASSFacility(log, creationtime, lifetime, 4) {
  DBGL;
  (*this).SetName("R_Reactor.");

  fIsStarted = false;
  fIsShutDown = false;
  fIsAtEndOfCycle = false;

  fStorage = 0;
  fFabricationPlant = 0;

  fFixedFuel = true;
  fIsStorage = false;

  fOutBackEndFacility = Pool;

  fCapacityFactor = CapacityFactor;
  fPower = power * fCapacityFactor;
  fEfficiencyFactor = 0.33;
  fElectricPower = fEfficiencyFactor * fPower;

  fHeavyMetalMass = HMMass;

  fBurnUp = -1;
  fCycleTime = (-1);

  fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
  fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
  fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
      (cSecond)(fCycleTime / fEvolutionDB.GetPower() * fPower));

  fReactorScheduler = 0;

  INFO << " A Reactor has been define :" << endl;
  INFO << "\t Fuel Composition is fixed (for now)! " << endl;
  INFO << "\t Creation time set at \t "
       << ((double)GetCreationTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t Life time (Operating's Duration) set at \t "
       << ((double)GetLifeTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t The Effective Thermal Power is \t " << (double)(fPower * 1e-6)
       << " MW (with Full Power " << power << " and " << CapacityFactor
       << " capacity factor)" << endl;
  INFO << "\t The Heavy Metal Mass in the Core set at "
       << (double)(fHeavyMetalMass) << " tons" << endl
       << endl;

  DBGL;
}
//________________________________________________________________________
Reactor::Reactor(CLASSLogger* log, FabricationPlant* fabricationplant,
                 CLASSBackEnd* Pool, cSecond creationtime, cSecond lifetime,
                 double Power, double HMMass, double CapacityFactor)
    : CLASSFacility(log, creationtime, lifetime, 4) {
  DBGL;
  (*this).SetName("R_Reactor.");

  fStorage = 0;
  fIsStarted = false;
  fIsShutDown = false;
  fIsAtEndOfCycle = false;

  fFabricationPlant = fabricationplant;
  fFixedFuel = false;

  fOutBackEndFacility = Pool;

  fBurnUp = -1;
  fHeavyMetalMass = HMMass;
  fCapacityFactor = CapacityFactor;
  fPower = Power * fCapacityFactor;
  fEfficiencyFactor = 0.33;
  fElectricPower = fEfficiencyFactor * fPower;

  fCycleTime = -1;  // BU in GWd/t

  fReactorScheduler = 0;

  INFO << " A Reactor has been define :" << endl;
  INFO << "\t Fuel Composition is not fixed (for now)! " << endl;
  INFO << "\t Creation time set at \t "
       << ((double)GetCreationTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t Life time (Operating's Duration) set at \t "
       << ((double)GetLifeTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t The Effective Thermal Power is \t " << (double)(fPower * 1e-6)
       << " MW (with Full Power " << Power << " and " << CapacityFactor
       << " capacity factor)" << endl;
  INFO << "\t The Heavy Metal Mass in the Core set at "
       << (double)(fHeavyMetalMass) << " tons" << endl
       << endl;

  DBGL;
}
//________________________________________________________________________
Reactor::Reactor(CLASSLogger* log, PhysicsModels* fueltypeDB,
                 FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
                 cSecond creationtime, cSecond lifetime, double Power,
                 double HMMass, double BurnUp, double CapacityFactor)
    : CLASSFacility(log, creationtime, lifetime, 4) {
  DBGL;
  (*this).SetName("R_Reactor.");

  fStorage = 0;
  fIsStarted = false;
  fIsShutDown = false;
  fIsAtEndOfCycle = false;

  fFabricationPlant = fabricationplant;

  fFixedFuel = false;

  fOutBackEndFacility = Pool;

  fBurnUp = BurnUp;
  fHeavyMetalMass = HMMass;
  fCapacityFactor = CapacityFactor;
  fPower = Power * fCapacityFactor;
  fEfficiencyFactor = 0.33;
  fElectricPower = fEfficiencyFactor * fPower;
  fCycleTime = (cSecond)(fBurnUp * 1e9 / (fPower)*fHeavyMetalMass * 3600 *
                         24);  // BU in GWd/t

  fReactorScheduler = new ReactorScheduler(log);
  fReactorScheduler->AddEntry(creationtime, new ReactorModel(fueltypeDB),
                              fBurnUp, fPower, fHeavyMetalMass);
  fSchedulePowerEvolution.insert(pair<cSecond, double>(creationtime, fPower));
  fScheduleHMMassEvolution.insert(
      pair<cSecond, double>(creationtime, fHeavyMetalMass));

  CheckListConsistency(fueltypeDB, fabricationplant);

  INFO << " A Reactor has been define :" << endl;
  INFO << "\t Fuel Composition is not fixed ! " << endl;
  INFO << "\t Creation time set at \t " << (double)(GetCreationTime() / cYear)
       << " year" << endl;
  INFO << "\t Life time (Operating's Duration) set at \t "
       << ((double)GetLifeTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t The Effective Thermal Power is \t " << (double)(fPower * 1e-6)
       << " MW (with Full Power " << Power << " and " << CapacityFactor
       << " capacity factor)" << endl;
  INFO << "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp)
       << " GWj/t" << endl;
  INFO << "\t The corresponding Cycle Time is\t "
       << ((double)fCycleTime) / ((double)cYear) << " year" << endl;
  INFO << "\t The Heavy Metal Mass in the Core set at "
       << (double)(fHeavyMetalMass) << " tons" << endl
       << endl;

  DBGL;
}

Reactor::Reactor(CLASSLogger* log, PhysicsModels* fueltypeDB,
                 FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
                 cSecond creationtime, cSecond lifetime, cSecond cycletime,
                 double HMMass, double BurnUp)
    : CLASSFacility(log, creationtime, lifetime, cycletime, 4) {
  DBGL;
  (*this).SetName("R_Reactor.");

  fIsStarted = false;
  fIsShutDown = false;
  fIsAtEndOfCycle = false;

  fStorage = 0;

  fFabricationPlant = fabricationplant;
  fFixedFuel = false;
  fBurnUp = BurnUp;
  fHeavyMetalMass = HMMass;

  fOutBackEndFacility = Pool;
  fPower = BurnUp * 3600. * 24. / (fCycleTime)*HMMass * 1e9;  // BU in GWd/t
  fEfficiencyFactor = 0.33;
  fElectricPower = fEfficiencyFactor * fPower;

  fReactorScheduler = new ReactorScheduler(log);
  fReactorScheduler->AddEntry(creationtime, new ReactorModel(fueltypeDB),
                              fBurnUp, fPower, fHeavyMetalMass);
  fSchedulePowerEvolution.insert(pair<cSecond, double>(creationtime, fPower));
  fScheduleHMMassEvolution.insert(
      pair<cSecond, double>(creationtime, fHeavyMetalMass));

  CheckListConsistency(fueltypeDB, fabricationplant);

  INFO << " A Reactor has been define :" << endl;
  INFO << "\t Fuel Composition is not fixed ! " << endl;
  INFO << "\t Creation time set at \t "
       << ((double)GetCreationTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t Life time (Operating's Duration) set at \t "
       << ((double)GetCreationTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t The Cycle Time set at\t "
       << ((double)fCycleTime) / ((double)cYear) << " year" << endl;
  INFO << "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp)
       << " GWj/t" << endl;
  INFO << "\t The corresponding Effective Thermal Power is \t "
       << (double)(fPower * 1e-6) << " MW" << endl;
  INFO << "\t The Heavy Metal Mass in the Core set at "
       << (double)(fHeavyMetalMass) << " tons" << endl
       << endl;

  DBGL;
}

Reactor::Reactor(CLASSLogger* log, EvolutionData* evolutivedb,
                 CLASSBackEnd* Pool, cSecond creationtime, cSecond lifetime,
                 double power, double HMMass, double BurnUp,
                 double CapacityFactor)
    : CLASSFacility(log, creationtime, lifetime, 4) {
  DBGL;
  (*this).SetName("R_Reactor.");

  fIsStarted = false;
  fIsShutDown = false;
  fIsAtEndOfCycle = false;

  fStorage = 0;
  fFabricationPlant = 0;

  fFixedFuel = true;
  fIsStorage = false;

  fOutBackEndFacility = Pool;

  fCapacityFactor = CapacityFactor;
  fPower = power * fCapacityFactor;
  fEfficiencyFactor = 0.33;
  fElectricPower = fEfficiencyFactor * fPower;

  fHeavyMetalMass = HMMass;

  double M0 = cZAIMass.GetMass(
      evolutivedb->GetIsotopicVectorAt(0.).GetActinidesComposition());

  fEvolutionDB = (*evolutivedb) * (fHeavyMetalMass / M0);

  fBurnUp = BurnUp;
  fCycleTime = (cSecond)(fBurnUp * 1e9 / (fPower)*fHeavyMetalMass * 3600 * 24);

  fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
  fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
  fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
      (cSecond)(fCycleTime / fEvolutionDB.GetPower() * fPower));

  fReactorScheduler = new ReactorScheduler(log);
  fReactorScheduler->AddEntry(creationtime, new ReactorModel(evolutivedb),
                              fBurnUp, fPower, fHeavyMetalMass);
  fSchedulePowerEvolution.insert(pair<cSecond, double>(creationtime, fPower));
  fScheduleHMMassEvolution.insert(
      pair<cSecond, double>(creationtime, fHeavyMetalMass));

  INFO << " A Reactor has been define :" << endl;
  INFO << "\t Fuel Composition is fixed ! " << endl;
  INFO << "\t Creation time set at \t " << (double)(GetCreationTime() / cYear)
       << " year" << endl;
  INFO << "\t Life time (Operating's Duration) set at \t "
       << ((double)GetLifeTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t The Cycle Time set at\t "
       << ((double)fCycleTime) / ((double)cYear) << " year" << endl;
  INFO << "\t The Effective Thermal Power is \t " << (double)(fPower * 1e-6)
       << " MW (with Full Power " << power << " and " << CapacityFactor
       << " capacity factor)" << endl;
  INFO << "\t The Heavy Metal Mass in the Core set at "
       << (double)(fHeavyMetalMass) << " tons" << endl
       << endl;

  DBGL;
}

Reactor::Reactor(CLASSLogger* log, EvolutionData* evolutivedb,
                 CLASSBackEnd* Pool, cSecond creationtime, cSecond lifetime,
                 cSecond cycletime, double HMMass, double BurnUp)
    : CLASSFacility(log, creationtime, lifetime, cycletime, 4) {
  DBGL;
  (*this).SetName("R_Reactor.");

  fIsStarted = false;
  fIsShutDown = false;
  fIsAtEndOfCycle = false;

  fStorage = 0;
  fFabricationPlant = 0;

  fFixedFuel = true;
  fIsStorage = false;

  fOutBackEndFacility = Pool;

  fPower = BurnUp * 3600. * 24. / (fCycleTime)*HMMass * 1e9;  // BU in GWd/t
  fCapacityFactor = 1;
  fEfficiencyFactor = 0.33;
  fElectricPower = fEfficiencyFactor * fPower;

  fHeavyMetalMass = HMMass;

  double M0 = cZAIMass.GetMass(
      evolutivedb->GetIsotopicVectorAt(0.).GetActinidesComposition());

  fEvolutionDB = (*evolutivedb) * (fHeavyMetalMass / M0);

  fBurnUp = BurnUp;

  fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
  fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
  fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
      (cSecond)(fCycleTime / fEvolutionDB.GetPower() * fPower));

  fReactorScheduler = new ReactorScheduler(log);
  fReactorScheduler->AddEntry(creationtime, new ReactorModel(evolutivedb),
                              fBurnUp, fPower, fHeavyMetalMass);

  fSchedulePowerEvolution.insert(pair<cSecond, double>(creationtime, fPower));
  fScheduleHMMassEvolution.insert(
      pair<cSecond, double>(creationtime, fHeavyMetalMass));

  INFO << " A Reactor has been define :" << endl;
  INFO << "\t Fuel Composition is fixed ! " << endl;
  INFO << "\t Creation time set at \t "
       << ((double)GetCreationTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t Life time (Operating's Duration) set at \t "
       << ((double)GetLifeTime()) / ((double)cYear) << " year" << endl;
  INFO << "\t The Cycle Time set at\t "
       << ((double)fCycleTime) / ((double)cYear) << " year" << endl;
  INFO << "\t The Effective Thermal Power is \t " << (double)(fPower * 1e-6)
       << " MW (with Full Power " << fPower << endl;
  INFO << "\t The Heavy Metal Mass in the Core set at "
       << (double)(fHeavyMetalMass) << " tons" << endl
       << endl;

  DBGL;
}

//________________________________________________________________________
Reactor::~Reactor() {}

//________________________________________________________________________
void Reactor::SetCycleTime(double cycletime) {
  fCycleTime = (cSecond)cycletime;

  if (fFixedFuel == true) {
    fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
        fCycleTime / fEvolutionDB.GetPower() * fPower);
    fBurnUp = fPower * fCycleTime / 3600. / 24. / fHeavyMetalMass;
  } else {
    fBurnUp = fPower * fCycleTime / 3600. / 24. / fHeavyMetalMass;
  }
}
//________________________________________________________________________
void Reactor::SetPower(double Power) {
  fPower = Power;

  if (fFixedFuel == true) {
    fCycleTime =
        (cSecond)(fBurnUp * 1e9 / (fPower)*fHeavyMetalMass * 3600 * 24);
    fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
        (cSecond)(fCycleTime / fEvolutionDB.GetPower() * fPower));
  } else
    fCycleTime = (cSecond)(fBurnUp * 1e9 / (fPower)*fHeavyMetalMass * 3600 *
                           24);  // BU in GWd/t
}
//________________________________________________________________________
void Reactor::SetBurnUp(double BU) {
  fBurnUp = BU;

  if (fFixedFuel == true) {
    fCycleTime =
        (cSecond)(fBurnUp * 1e9 / (fPower)*fHeavyMetalMass * 3600 * 24);
    fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
        (cSecond)(fCycleTime / fEvolutionDB.GetPower() * fPower));
  } else
    fCycleTime =
        (cSecond)(fBurnUp * 1e9 / (fPower)*fHeavyMetalMass * 3600 * 24);
}

//________________________________________________________________________
void Reactor::SetEvolutionDB(EvolutionData evolutionDB) {
  DBGL;

  // fEvolutionDB.DeleteEvolutionDataCopy();

  double M0 = cZAIMass.GetMass(
      evolutionDB.GetIsotopicVectorAt(0.).GetActinidesComposition());
  fEvolutionDB = evolutionDB * (fHeavyMetalMass / M0);

  fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(
      (cSecond)(fCycleTime / fEvolutionDB.GetPower() * fPower));
  fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);

  DBGL;
}

//________________________________________________________________________
void Reactor::SetNewFuel(EvolutionData ivdb) {
  DBGL;
  SetEvolutionDB(ivdb);
  DBGL;
}

//________________________________________________________________________
void Reactor::Evolution(cSecond t) {
  DBGL;

  if (fIsShutDown || t < GetCreationTime())
    return;  // Reactor stop or not started...

  if (Norme(fInsideIV) != 0 && fIsStarted) {
#pragma omp critical(ParcPowerUpdate)
    { GetParc()->AddToPower(fPower, fElectricPower); }
  } else if (fIsStarted) {
    WARNING << " Reactor should be working but have no Heavy Nucleus Inside. "
               "It's not working so have a zero power..."
            << " Time : " << t / cYear << " years" << endl;
  }

  if (t < fInternalTime) return;
  if (t == fInternalTime && t != GetCreationTime()) return;

  if (t == GetCreationTime() && !fIsStarted)  // Start of the Reactor
  {
    fIsAtEndOfCycle = true;
    fInsideIV = fIVBeginCycle;
    fInternalTime = t;
    fInCycleTime = 0;
  }

  // Check if the Reactor if started ...
  if (!fIsStarted)
    return;  // If the reactor just start don't need to make Fuel evolution

  cSecond EvolutionTime =
      t - fInternalTime;  // Calculation of the evolution time (relativ)

  if (EvolutionTime + fInCycleTime == fCycleTime)  // End of Cycle
  {
    fIsAtEndOfCycle = true;
    fInternalTime += EvolutionTime;  // Update Internal Time
    fInCycleTime += EvolutionTime;   // Update InCycleTime

    fInsideIV = fEvolutionDB.GetIsotopicVectorAt(
        (cSecond)(fInCycleTime / fEvolutionDB.GetPower() *
                  fPower));  // update the fuel composition
    if (t >= GetCreationTime() + GetLifeTime())
      fIsShutDown = true;  // if the Next Cycle don't 'Exist...

  } else if (EvolutionTime + fInCycleTime < fCycleTime)  // During Cycle
  {
    fInternalTime += EvolutionTime;  // Update Internal Time
    fInCycleTime += EvolutionTime;   // Update InCycleTime

    fInsideIV = fEvolutionDB.GetIsotopicVectorAt(
        (cSecond)(fInCycleTime / fEvolutionDB.GetPower() *
                  fPower));  // update the fuel composition
    if (t >= GetCreationTime() + GetLifeTime()) fIsShutDown = true;
  } else {
    // Evolution goes after the end of cycle.... check it
    ERROR << " " << (*this).GetName() << endl;
    ERROR << " Evolution Time is " << t << endl;
    ERROR
        << " Evolution is too long! There is a problem in Reactor evolution at "
        << t / cYear << endl;
    ERROR << " This is too long of : "
          << EvolutionTime + fInCycleTime - fCycleTime << endl;
    ERROR << " I have spend " << fInCycleTime + EvolutionTime
          << " and should have been " << fCycleTime << endl;
    exit(1);
  }

  DBGL;
}

//________________________________________________________________________
void Reactor::Dump() {
  DBGL;

  if (fInternalTime < GetCreationTime()) return;
  if (fIsShutDown && !fIsStarted) return;  // Reactor stopped...
  if (!fIsAtEndOfCycle && !fIsShutDown) return;

  // First trash the irradiated fuel
  if (fIsAtEndOfCycle && !fIsShutDown) {
    if (fIsStarted)  // A Cycle has already been done
    {
      fOutBackEndFacility->AddIV(fInsideIV);
      AddCumulativeIVOut(fInsideIV);
    } else
      fIsStarted = true;  // Just start the first cycle

  } else if (fIsAtEndOfCycle && fIsShutDown)  // shutdown at end of Cycle
  {
    fOutBackEndFacility->AddIV(fIVOutCycle);
    AddCumulativeIVOut(fIVOutCycle);
    fInsideIV.Clear();
    fInCycleTime = 0;
    fIsStarted = false;                        // shut down the Reactor
  } else if (!fIsAtEndOfCycle && fIsShutDown)  // shutdown during Cycle
  {
    fOutBackEndFacility->AddIV(fInsideIV);
    AddCumulativeIVOut(fInsideIV);
    fInsideIV.Clear();
    fInCycleTime = 0;
    fIsStarted = false;  // shut down the Reactor
  }

  DBGL;

  // Get the new Fuel !
  ScheduleEntry* NextFuel = fReactorScheduler->GetEntryAt(fInternalTime);
  fPower = NextFuel->GetPower();
  fHeavyMetalMass = NextFuel->GetHeavyMetalMass();
  fElectricPower = fPower * fEfficiencyFactor;
  SetBurnUp(NextFuel->GetBurnUp());

  fPowerEvolution.insert(pair<cSecond, double>(fInternalTime, fPower));
  fHMMassEvolution.insert(
      pair<cSecond, double>(fInternalTime, fHeavyMetalMass));

  if (NextFuel->GetReactorModel()->GetPhysicsModels())
    fFixedFuel = false;
  else if (NextFuel->GetReactorModel()->GetEvolutionData())
    fFixedFuel = true;
  else {
    ERROR << typeid(NextFuel->GetReactorModel()).name() << endl;
    ERROR << "WRONG Fuel Format Correct it !! " << endl;
    exit(1);
  }
  DBGL;

  if (fFixedFuel) {
    DBGL;
    if (fIsAtEndOfCycle && !fIsShutDown) {
      SetEvolutionDB(*(NextFuel->GetReactorModel()->GetEvolutionData()));
      fIsAtEndOfCycle = false;

      if (!GetParc()->GetStockManagement() && fIsStorage) {
        IsotopicVector BuildIVtmp;
        IsotopicVector OutIncomePart;

        // Get The Storage Compostion
        BuildIVtmp.Add(fStorage->GetInsideIV().GetIsotopicQuantity());
        // Get the rest after IVIn creation
        BuildIVtmp -= fIVInCycle;
        // Get the OutIncome part form this rest
        OutIncomePart.Add(BuildIVtmp.GetIsotopicQuantityNeeded());
        // Take what you can from Storage...
        fStorage->TakeFromStock(fIVInCycle - OutIncomePart);
        // And Get the rest from OutIncome
        GetParc()->AddOutIncome(OutIncomePart);

      } else
        GetParc()->AddOutIncome(fIVInCycle);

      fInsideIV = fIVBeginCycle;
      AddCumulativeIVIn(fIVBeginCycle);

      fInCycleTime = 0;
    }
    DBGL;
  } else {
    DBGL;
    if (!GetParc()->GetStockManagement()) {
      ERROR << " Can't have unfixedFuel without stock management" << endl;
      ERROR << " Add in your input : gCLASS->SetStockManagement(true);" << endl;
      exit(1);
    }

    if (fIsAtEndOfCycle && !fIsShutDown) {
      fIsAtEndOfCycle = false;
      SetNewFuel(fFabricationPlant->GetReactorEvolutionDB(GetId()));
      fFabricationPlant->TakeReactorFuel(GetId());

      fInsideIV = fIVBeginCycle;
      AddCumulativeIVIn(fIVBeginCycle);

      fInCycleTime = 0;
    }

    DBGL;
  }

  DBGL;
}

//________________________________________________________________________
cSecond Reactor::GetNextCycleTime(cSecond time) {
  DBGL;
  cSecond LastCycle = fInternalTime - fInCycleTime;

  while (LastCycle < time) {
    ScheduleEntry* entry = fReactorScheduler->GetEntryAt(LastCycle);
    double BU = entry->GetBurnUp();
    double Power = entry->GetPower();
    double HMMass = entry->GetHeavyMetalMass();

    cSecond cycletime = (cSecond)(BU * 1e9 / (Power)*HMMass * 3600 * 24);
    LastCycle += cycletime;
  }

  DBGL;
  return LastCycle;
}

//________________________________________________________________________
void Reactor::CheckListConsistency(PhysicsModels* fueltypeDB,
                                   FabricationPlant* fabricationplant) {
  // Get Lists containers defined in EqM and FP
  map<string, IsotopicVector> StreamList =
      fueltypeDB->GetEQM()->GetAllStreamList();
  map<string, vector<Storage*> > Stocks = fabricationplant->GetAllStorage();

  vector<string> StreamListEqM;
  vector<string> StreamListFP;

  vector<bool> CheckOnLists;
  bool AreListsConsistent = true;

  // Iterators
  map<string, vector<Storage*> >::iterator it_s_vS;
  map<string, IsotopicVector>::iterator it_s_IV;

  for (it_s_vS = Stocks.begin(); it_s_vS != Stocks.end(); it_s_vS++)
    StreamListFP.push_back((*it_s_vS).first);

  for (it_s_IV = StreamList.begin(); it_s_IV != StreamList.end(); it_s_IV++)
    StreamListEqM.push_back((*it_s_IV).first);

  if ((int)StreamListEqM.size() != (int)StreamListFP.size()) {
    WARNING << " Not the same number of lists in FP and associated EqM" << endl;
    AreListsConsistent = false;
  }

  else {
    for (int i = 0; i < (int)StreamListEqM.size(); i++) {
      CheckOnLists.push_back(false);
    }

    for (int i = 0; i < (int)StreamListEqM.size(); i++) {
      for (int j = 0; j < (int)StreamListFP.size(); j++) {
        if (StreamListEqM[i] == StreamListFP[j]) {
          CheckOnLists[i] = true;
        }
      }
    }

    for (int i = 0; i < (int)CheckOnLists.size(); i++) {
      if (!CheckOnLists[i]) {
        AreListsConsistent = false;
      }
    }
  }

  if (!AreListsConsistent) {
    ERROR << "Lists defined in Fabrication plant and in EqM are not the same."
          << endl;
    ERROR << "Lists defined in Fabrication plant :" << endl;
    for (int j = 0; j < (int)StreamListFP.size(); j++) {
      ERROR << "-> " << StreamListFP[j] << "  " << endl;
    }
    ERROR << "Lists defined in EqM :" << endl;

    for (int i = 0; i < (int)StreamListEqM.size(); i++) {
      ERROR << "-> " << StreamListEqM[i] << "  " << endl;
    }
    ERROR << "Check in your scenario." << endl;

    exit(1);
  }
}
