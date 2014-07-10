#include "PhysicModels.hxx"

//________________________________________________________________________
//
//		PhysicModels
//
//
//
//
//________________________________________________________________________
PhysicModels::PhysicModels(XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM ):CLASSObject()
{

		fXSModel		= XS;
		fEquivalenceModel	= EM;
		fIrradiationModel	= IM;


}
//________________________________________________________________________
PhysicModels::PhysicModels(CLASSLogger* log, XSModel* XS, EquivalenceModel* EM, IrradiationModel* IM ):CLASSObject(log)
{

	fXSModel		= XS;
	fEquivalenceModel	= EM;
	fIrradiationModel	= IM;


}

//________________________________________________________________________
EvolutionData PhysicModels::GenerateEvolutionData(IsotopicVector IV, double cycletime, double Power)
{

	return fIrradiationModel->GenerateEvolutionData(IV, fXSModel->GetCrossSections(IV), Power, cycletime);
}
//________________________________________________________________________
