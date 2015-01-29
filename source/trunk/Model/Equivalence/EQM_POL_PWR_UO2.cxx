#include "EquivalenceModel.hxx"
#include "EQM_POL_PWR_UO2.hxx"
#include "CLASSLogger.hxx"


// ________________________________________________________________________
// EQM_POL_PWR_UO2
//
// Brief description
// ________________________________________________________________________




//Constructor(s)
EQM_POL_PWR_UO2::EQM_POL_PWR_UO2()
{
	
	// Fertile
	ZAI U8(92 ,238 ,0) ;
	fFertileList = U8*1;
	// Fissile
	ZAI U5(92 ,235 ,0) ;
	// ...
	fFissileList = U5*1;
}
// _______________________________________________________________________
double EQM_POL_PWR_UO2::GetFissileMolarFraction ( IsotopicVector Fissil , IsotopicVector Fertil , double BurnUp )
{

	return 0.0125575 + 0.00053602*BurnUp + 1.80023e-06*BurnUp*BurnUp ;
	
}