#ifndef _EQM_LIN_PWR_MOX_HXX
#define _EQM_LIN_PWR_MOX_HXX

#include "EquivalenceModel.hxx"

#include <string>

/*!
 \file
 \brief Header file for EQM_LIN_PWR_MOX class.


 @author BaM
 @version 1.0
 */

using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a EQM_LIN_PWR_MOX.
 The aim of these class is to constuct a fuel from an equivalence model
 based on a Linear Eq Model @f$BU = \alpha_{0} + \sum_{i\in fissile}alpha_{i}\cdot n_{i} @f$
 For one set of  @f$\alpha @f$ values the fabrication time is fixed in order to decrease 
 the dimentionality by 1 (@f$^{241}Am@f$ is thus fixed by this decay fixed decay time) .
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_LIN_PWR_MOX : public EquivalenceModel
{
	public :
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// Simple constructor
	
	/*!
	 Make a new EQM_LIN_PWR_MOX
	 \param WeightPath : Path to the file containing the @f$\alpha_{i}@f$. The file is format as :
	 
	 @f$\alpha_{0}\,\alpha_{^{238}Pu}\,\alpha_{^{239}Pu}\,\alpha_{^{240}Pu}\,\alpha_{^{241}Pu}\alpha_{^{242}Pu}@f$
	 
	 */
	EQM_LIN_PWR_MOX(string WeightPath);
	//}
	
	//{
	/// Logger constructor
	
	/*!
	 Make a new EQM_LIN_PWR_MOX
	 \param log : use for the log
	 \param WeightPath : Path to the file containing the @f$\alpha_{i}@f$. The file is format as :
	 
	 @f$\alpha_{0}\,\alpha_{^{238}Pu}\,\alpha_{^{239}Pu}\,\alpha_{^{240}Pu}\,\alpha_{^{241}Pu}\alpha_{^{242}Pu}@f$
	 
	 */
	EQM_LIN_PWR_MOX(CLASSLogger* log, string WeightPath);
	
	~EQM_LIN_PWR_MOX();

	virtual vector<double> BuildFuel(double BurnUp, double HMMass, vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray );

	private :

	string fWeightPath;		//!< The full path to the file containing the @f$\alpha_{i}@f$
	vector<double> fFuelParameter;	//!< The vector of @f$\alpha_{i}@f$

};

#endif

