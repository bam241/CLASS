#ifndef _EQM_PWR_POL_UO2_HXX
#define _EQM_PWR_POL_UO2_HXX

#include "EquivalenceModel.hxx"

using namespace std;

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−--------------−−−−−−−−−−−−−−−//
//! Define an EquivalenceModel based on a polynomial fit

/*!
 Defines a EQM_PWR_POL_UO2
 It returns the @f$^{235}U@f$ enrichment e according to this polynom :
 
 @f$e = \alpha_{0} + \alpha_{1}\cdot Burnup + \alpha_{2}\cdot Burnup \cdot Burnup @f$
 
 BU : Maximum achievable burnup
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EQM_PWR_POL_UO2 : public EquivalenceModel
{
	
public:
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// normal constructor
	/*!
	 Create a EQM_PWR_POL_UO2
	 \param  PathToWeightFile :  Path to the file containing the  @f$\alpha_{i}@f$
	 Format :  @f$PARAM@f$  @f$\alpha_{0}@f$  @f$\alpha_{1}@f$  @f$\alpha_{2}@f$
	 */
	EQM_PWR_POL_UO2(string PathToWeightFile);
	//}
	
	//{
	/// logger constructor
	/*!
	 Create a EQM_PWR_POL_UO2
	 \param log : Use for the log
	 \param  PathToWeightFile :  Path to the file containing the  @f$\alpha_{i}@f$
	 Format :  @f$PARAM@f$  @f$\alpha_{0}@f$  @f$\alpha_{1}@f$  @f$\alpha_{2}@f$
	 */
	EQM_PWR_POL_UO2(CLASSLogger* log, string PathToWeightFile);
	//}
	
	//@}
	/**This function IS the equivalence model**/
	map < string , double> GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp);
		
	
private:
	
	void ReadWeightFile(string PathToWeightFile); //!< Function to read the weight file & file the parameters
	
	double fParam_Bu_0 ;		//!<  @f$\alpha_{0}@f$
	double fParam_Bu ;		//!<  @f$\alpha_{1}@f$
	double fParam_BuSquare ;	//!<  @f$\alpha_{2}@f$

};

#endif