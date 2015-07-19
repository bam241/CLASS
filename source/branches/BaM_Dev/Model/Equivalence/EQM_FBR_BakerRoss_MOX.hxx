#ifndef _EQM_FBR_BakerRoss_MOX_HXX
#define _EQM_FBR_BakerRoss_MOX_HXX

#include "EquivalenceModel.hxx"

#include <string>

using namespace std;

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on Baker and Ross formula

/*!
 The aim of these class is to constuct a fuel from an equivalence model
 based on a @f$^{239}Pu@f$ equivalent from Baker and Ross formula :
 
 It returns Pu content (E) needed for the FBR-Na loaded with a given Pu vector according
  :
 
 @f$E = \frac{E_{ref} - \sum_{fertile}N_{i}W_{i} }{\sum_{fissile}N_{i}W_{i}-\sum_{fertile}N_{i}W_{i}}@f$
 
 with :
 
 @f$W_i = \frac{\alpha_{i} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
 and
 @f$\alpha_{i} = \bar{\nu_{i}}\cdot\sigma_{i}^{fis} - \sigma_{i}^{abs}@f$

 @author BLG
 @version 3.0
 */
//________________________________________________________________________

class EQM_FBR_BakerRoss_MOX : public EquivalenceModel
{
	public :
	
	/*!
	 \name Constructor
	 */
	//@{
	//{
	/// Logger constructor
	
	/*!
	 Make a new EQM_FBR_BakerRoss_MOX : the default values have been calculated for a FBR-Na of type : ESFR like (without blanket)
	 \param log : use for the log
	 \param Weight_U_235  : reactivity weight @f$W_{^{235}U} = \frac{\alpha_{{^{235}U}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_238 : reactivity weight @f$W_{^{238}Pu} = \frac{\alpha_{{^{238}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_240 : reactivity weight @f$W_{^{240}Pu} = \frac{\alpha_{{^{240}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_241 : reactivity weight @f$W_{^{241}Pu} = \frac{\alpha_{{^{241}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_242 : reactivity weight @f$W_{^{242}Pu} = \frac{\alpha_{{^{242}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Am_241 : reactivity weight @f$W_{^{241}Pu} = \frac{\alpha_{{^{241}Am}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param EquivalentFissile : reference fresh fuel @f$^{239}Pu@f$ content neeed in a @f$^{238}U@f$ + @f$^{239}Pu@f$ fuel to satisfy criticality at t = 0
	 */
	EQM_FBR_BakerRoss_MOX(CLASSLogger* log, double Weight_U_235 = 0.791135, double Weight_Pu_238 = 0.686385, double Weight_Pu_240 = 0.13553, double Weight_Pu_241 = 1.54572 , double Weight_Pu_242 = 0.0829001, double Weight_Am_241 = -0.336945, double EquivalentFissile = 0.103213);
	//}
	
	//{
	/// normal constructor
	
	/*!
	 Make a new EQM_FBR_BakerRoss_MOX : the default values have been calculated for a FBR-Na of type : ESFR like (without blanket)
	 \param Weight_U_235  : reactivity weight @f$W_{^{235}U} = \frac{\alpha_{{^{235}U}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_238 : reactivity weight @f$W_{^{238}Pu} = \frac{\alpha_{{^{238}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_240 : reactivity weight @f$W_{^{240}Pu} = \frac{\alpha_{{^{240}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_241 : reactivity weight @f$W_{^{241}Pu} = \frac{\alpha_{{^{241}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Pu_242 : reactivity weight @f$W_{^{242}Pu} = \frac{\alpha_{{^{242}Pu}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param Weight_Am_241 : reactivity weight @f$W_{^{241}Pu} = \frac{\alpha_{{^{241}Am}} - \alpha_{^{238}U} }{\alpha_{^{239}Pu}-\alpha_{^{238}U}}@f$
	 \param EquivalentFissile : reference fresh fuel @f$^{239}Pu@f$ content neeed in a @f$^{238}U@f$ + @f$^{239}Pu@f$ fuel to satisfy criticality at t = 0
	 */
	EQM_FBR_BakerRoss_MOX(double Weight_U_235 = 0.791135, double Weight_Pu_238 = 0.686385, double Weight_Pu_240 = 0.13553, double Weight_Pu_241 = 1.54572 , double Weight_Pu_242 = 0.0829001, double Weight_Am_241 = -0.336945, double EquivalentFissile = 0.103213);
	//}
	//@}
	
	virtual double GetFissileMolarFraction(IsotopicVector Fissil, IsotopicVector Fertil, double BurnUp = 0);

	private :
	double fReferenceFissilContent;   //!<The reference fraction of Pu for BurnUp = 100Gwj/t with a ideal model(only @f$^{239}Pu@f$ and &@f$^{238}U@f$ are taken into account)

	/*!
	 \name Reactivity coefficients :
	 */
	//@{
	double fWeight_U_235;	//!< weight for @f$^{235}U@f$
	double fWeight_Pu_238;	//!< weight for @f$^{238}Pu@f$
	double fWeight_Pu_240;	//!< weight for @f$^{240}Pu@f$
	double fWeight_Pu_241;	//!< weight for @f$^{241}Pu@f$
	double fWeight_Pu_242;	//!< weight for @f$^{242}Pu@f$
	double fWeight_Am_241;	//!< weight for @f$^{241}Am@f$
	//@}
};

#endif


