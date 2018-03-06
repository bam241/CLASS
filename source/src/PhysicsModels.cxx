#include "PhysicsModels.hxx"

//________________________________________________________________________
//
//		PhysicsModels
//
//
//
//
//________________________________________________________________________

PhysicsModels::PhysicsModels() : CLASSObject() {
  fXSM = 0;
  fEQM = 0;
  fIM = 0;
}
//________________________________________________________________________
PhysicsModels::PhysicsModels(XSModel* XS, EquivalenceModel* EM,
                             IrradiationModel* IM)
    : CLASSObject() {
  fXSM = XS;
  fEQM = EM;
  fIM = IM;

  int Z_ZAIThreshold = fIM->GetZAIThreshold();
  fXSM->SetZAIThreshold(Z_ZAIThreshold);
}
//________________________________________________________________________
PhysicsModels::PhysicsModels(CLASSLogger* log, XSModel* XS,
                             EquivalenceModel* EM, IrradiationModel* IM)
    : CLASSObject(log) {
  fXSM = XS;
  fEQM = EM;
  fIM = IM;

  int Z_ZAIThreshold = fIM->GetZAIThreshold();
  fXSM->SetZAIThreshold(Z_ZAIThreshold);
}

//________________________________________________________________________
EvolutionData PhysicsModels::GenerateEvolutionData(IsotopicVector IV,
                                                   double cycletime,
                                                   double Power) {
  fXSM->isIVInDomain(IV);

  return fIM->GenerateEvolutionData(
      IV, fXSM->GetCrossSections(IV), Power, cycletime);
}
//________________________________________________________________________
