#ifndef _XSMODEL_HXX
#define _XSMODEL_HXX


/*!
 \file
 \brief Header file for XSMODEL class.
 
 
 @author BLG
 @version 1.0
 */
#include "EvolutionData.hxx"
#include "CLASSObject.hxx"


using namespace std;

class IsotopicVector;

//-----------------------------------------------------------------------------//
/*!
 Define a XS Interpolator 
This is the class for method related to XS prediction

see @../Model/XS/XSM_CLOSEST.hxx 
see @../Model/XS/XSM_MLP.hxx 
...

!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!
Never instantiate XSModel in your CLASS input but it's derivated class


 @author BLG
 @version 1.0
 */
//________________________________________________________________________


class XSModel : public CLASSObject
{

	public : 

	XSModel();
	XSModel(CLASSLogger* log);
	/// Pure virtual method called to estimates mean cross sections.
	/*!
	 Estimates the mean cross sections evolution according the fresh fuel composition
	 \param vector<IsotopicVector> IV, Fresh fuel composition
	 \param double t =0, can be used in XSM_CLOSEST to recalculate 
	 distance between IV and DataBase at each time step (deprecated)
	 */
	virtual  EvolutionData GetCrossSections(IsotopicVector IV,double t=0) = 0 ;
	
	/// Check either the IsotopicVector IV is in the validity domain of the models.
	/*!
	 return true if IV is in ValidityDomain
	 return false + a warning if IV is not in ValidityDomain
	 \param vector<IsotopicVector> IV, Fresh fuel composition
	 */
	virtual  bool isIVInDomain(IsotopicVector IV) ;
	
	void SetZAIThreshold(int Z_Threshold){fZAIThreshold = Z_Threshold;}//!< Set the Z threshold : ZAI with Z < fZAIThreshold are not manage by CLASS
	int  GetZAIThreshold(){return fZAIThreshold;}//!< Get the Z threshold

	protected :
	map< ZAI, pair<double,double> > fFreshFuelDomain; //!< Fresh fuel range : map<ZAI<min edge ,max edge >>
	int fZAIThreshold;	//!< Z threshold for handling nuclei mean cross section (take only ZAI reaction of Z>=fZAIThresold)
};

#endif

