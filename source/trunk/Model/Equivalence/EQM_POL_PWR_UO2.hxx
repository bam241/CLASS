#ifndef _EQM_POL_PWR_UO2_HXX
#define _EQM_POL_PWR_UO2_HXX
#include "EquivalenceModel.hxx" 
using namespace std;
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−//
/*!
Define a EQM_POL_PWR_UO2
Explain briefly what is it . @author YourName
@version 3.0
*/
//________________________________________________________________________
class EQM_POL_PWR_UO2 : public EquivalenceModel
{
  public:
/*Constructor*/
EQM_POL_PWR_UO2(/*parameters*/ ); //!< Explain what is the parameters (if any)

/**This function IS the equivalence model**/
double GetFissileMolarFraction(IsotopicVector Fissil, IsotopicVector Fertil,double BurnUp) ; // !<Return the molar fraction of fissile element

  private:
/*Your private variables */ 

};

#endif