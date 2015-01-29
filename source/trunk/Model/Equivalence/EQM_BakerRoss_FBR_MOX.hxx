#ifndef _EQM_BakerRoss_FBR_MOX_HXX
#define _EQM_BakerRoss_FBR_MOX_HXX

#include "EquivalenceModel.hxx"

#include <string>

/*!
 \file
 \brief Header file for EQM_BakerRoss_FBR_MOX class.


 @author BLG
 @version 1.0
 */

using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a EQM_MLP_MOX.
 The aim of these class is to constuct a fuel from an equivalence model
 based on a 239Pu equivalent Model

 @author BLG
 @version 3.0
 */
//________________________________________________________________________

class EQM_BakerRoss_FBR_MOX : public EquivalenceModel
{
	public :
	EQM_BakerRoss_FBR_MOX(CLASSLogger* log, double Weight_U_235 =7.91642e-01  , double Weight_Pu_238=6.84741e-01  , double Weight_Pu_240=1.37467e-01  , double Weight_Pu_241=1.54546e+00  , double Weight_Pu_242=8.11544e-02  , double Weight_Am_241=-3.36839e-01  , double EquivalentFissile = 0.106);
	EQM_BakerRoss_FBR_MOX( double Weight_U_235 =7.91642e-01  , double Weight_Pu_238=6.84741e-01  , double Weight_Pu_240=1.37467e-01  , double Weight_Pu_241=1.54546e+00  , double Weight_Pu_242=8.11544e-02  , double Weight_Am_241=-3.36839e-01  , double EquivalentFissile = 0.106);

	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);

	private :
 	//The reference fraction of Pu for BurnUp=100Gwj/t with a ideal model(only pu239 and U238 are taken into account) 
	double fReferenceFissilContent;  

    // Reactivity coefficients for the isotopes
	double fWeight_U_235;
	double fWeight_Pu_238;
	double fWeight_Pu_240;
	double fWeight_Pu_241;
	double fWeight_Pu_242;
	double fWeight_Am_241;
};

#endif


