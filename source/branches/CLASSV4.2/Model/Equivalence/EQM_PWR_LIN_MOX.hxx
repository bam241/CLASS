#ifndef _EQM_PWR_LIN_MOX_HXX
#define _EQM_PWR_LIN_MOX_HXX

#include "EquivalenceModel.hxx"

#include <string>

using namespace std;

//-----------------------------------------------------------------------------//
//! Defines an EquivalenceModel based on a linear fit

/*!
 The aim of these class is to constuct a fuel from an equivalence model
 based on a Linear Eq Model @f$BU = \alpha_{0} + \sum_{i\in fissile}\alpha_{i}\cdot n_{i} @f$
 For one set of  @f$\alpha @f$ values the fabrication time is fixed in order to decrease 
 the dimentionality by 1 (@f$^{241}Am@f$ is thus fixed by this fixed decay time) .
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_PWR_LIN_MOX : public EquivalenceModel
{
	public :
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// Simple constructor
	
	/*!
	 Make a new EQM_PWR_LIN_MOX
	 \param WeightPath : Path to the file containing the @f$\alpha_{i}@f$. The file is format as :
	 
	 @f$PARAM@f$   @f$\alpha_{0}@f$  @f$\alpha_{^{238}Pu}@f$   @f$\alpha_{^{239}Pu}@f$  @f$\alpha_{^{240}Pu}@f$   @f$\alpha_{^{241}Pu}@f$   @f$\alpha_{^{242}Pu}@f$
	 
	 */
	EQM_PWR_LIN_MOX(string WeightPath);
	//}
	
	//{
	/// Logger constructor
	
	/*!
	 Make a new EQM_PWR_LIN_MOX
	 \param log : use for the log
	 \param WeightPath : Path to the file containing the @f$\alpha_{i}@f$. The file is format as :
	 
	  @f$PARAM@f$   @f$\alpha_{0}@f$  @f$\alpha_{^{238}Pu}@f$   @f$\alpha_{^{239}Pu}@f$  @f$\alpha_{^{240}Pu}@f$   @f$\alpha_{^{241}Pu}@f$   @f$\alpha_{^{242}Pu}@f$
	 
	 */
	EQM_PWR_LIN_MOX(CLASSLogger* log, string WeightPath);
	//}
	
	~EQM_PWR_LIN_MOX();
	//@}

	virtual vector<double> BuildFuel(double BurnUp, double HMMass, vector < vector <IsotopicVector> > StreamArray);

	private :

	string fWeightPath;		//!< The full path to the file containing the @f$\alpha_{i}@f$
	vector<double> fFuelParameter;	//!< The vector of @f$\alpha_{i}@f$

};

#endif

