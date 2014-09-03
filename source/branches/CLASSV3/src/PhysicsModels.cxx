#include "PhysicsModels.hxx"

//________________________________________________________________________
//
//		PhysicsModels
//
//
//
//
//________________________________________________________________________



PhysicsModels::PhysicsModels():CLASSFuel()
{

	fXSModel		= 0;
	fEquivalenceModel	= 0;
	fIrradiationModel	= 0;


}
//________________________________________________________________________
PhysicsModels::PhysicsModels(XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM ):CLASSFuel()
{

	fXSModel		= XS;
	fEquivalenceModel	= EM;
	fIrradiationModel	= IM;


}
//________________________________________________________________________
PhysicsModels::PhysicsModels(CLASSLogger* log, XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM ):CLASSFuel(log)
{

	fXSModel		= XS;
	fEquivalenceModel	= EM;
	fIrradiationModel	= IM;


}

//________________________________________________________________________
EvolutionData PhysicsModels::GenerateEvolutionData(IsotopicVector IV, double cycletime, double Power)
{

	return fIrradiationModel->GenerateEvolutionData(IV, fXSModel->GetCrossSections(IV), Power, cycletime);
}
//________________________________________________________________________
