#ifndef _EQM_QUAD_PWR_MOX_HXX
#define _EQM_QUAD_PWR_MOX_HXX

#include "EquivalenceModel.hxx"

#include <string>

/*!
 \file
 \brief Header file for EQM_MLP_MOX class.


 @author BaM
 @version 1.0
 */

using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a EQM_QUAD_PWR_MOX.
 The aim of these class is to constuct a fuel from an equivalence model
 based on a Quadratic Pu equivalent Model
 The Plutonium content e is calculated using :
 
 @f$ e = \alpha_{0} + \sum_{i\in Pu}^{N} \left(\alpha_{i} \cdot n_{i}\ + \sum_{j\leq i} \alpha_{ij} \cdot n_{i}\cdot n_{j}\right)@f$
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_QUAD_PWR_MOX : public EquivalenceModel
{
	public :

	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// normal constructor
	/*!
	 Create a EQM_POL_PWR_UO2
	 \param  WeightPath :  Path to the file containing the  @f$\alpha_{i}@f$
	 Format :  @f$PARAM@f$   @f$\alpha_{^{238}Pu}@f$  @f$\alpha_{^{238}Pu ^{238}Pu}@f$  @f$\alpha_{^{238}Pu ^{239}Pu}@f$ .... @f$\alpha_{^{238}Pu ^{242}Pu}@f$ @f$\alpha_{^{239}Pu}@f$ ...
	 @f$\alpha_{^{242}Pu}@f$  @f$\alpha_{^{242}Pu ^{242}Pu}@f$ @f$\alpha_{0}@f$
	 */
	EQM_QUAD_PWR_MOX(string WeightPath);
	//}
	
	//{
	/// Logger constructor
	/*!
	 Create a EQM_POL_PWR_UO2
	 \param log : Use for the log
	 \param  WeightPath :  Path to the file containing the  @f$\alpha_{i}@f$
	 Format :  @f$PARAM@f$   @f$\alpha_{^{238}Pu}@f$  @f$\alpha_{^{238}Pu ^{238}Pu}@f$  @f$\alpha_{^{238}Pu ^{239}Pu}@f$ .... @f$\alpha_{^{238}Pu ^{242}Pu}@f$ @f$\alpha_{^{239}Pu}@f$ ...
	 @f$\alpha_{^{242}Pu}@f$  @f$\alpha_{^{242}Pu ^{242}Pu}@f$ @f$\alpha_{0}@f$
	 */
	EQM_QUAD_PWR_MOX(CLASSLogger* log, string WeightPath);
	//}
	
	~EQM_QUAD_PWR_MOX();
	//@}
	
	virtual double GetFissileMolarFraction(IsotopicVector Fissil,IsotopicVector Fertil,double BurnUp);

	private :

	string fWeightPath;		//!< Path to the file containing the @f$\alpha_{ij}@f$
	vector<double> fFuelParameter;	//!< vector containing the @f$\alpha_{ij}@f$

};

#endif

