#ifndef _XSMODEL_HXX
#define _XSMODEL_HXX


/*!
 \file
 \brief Header file for XSMODEL class.
 
 
 @author BaM
 @author BLG
 @version 1.0
 */
#include "EvolutionData.hxx"
#include "CLASSObject.hxx"


using namespace std;

class IsotopicVector;

//-----------------------------------------------------------------------------//
//!  Defines a mean cross section predictor

/*!
This is the mother class for methods related to XS prediction

\warning
 Never instantiate XSModel in your CLASS input but it's derivated class

 @see XSM_CLOSEST
 @see XSM_MLP

 @author BaM
 @author BLG
 @version 1.0
 */
//________________________________________________________________________


class XSModel : public CLASSObject
{

	public : 

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	XSModel(); //!<Default constructor
	
	
	XSModel(CLASSLogger* log); //!<Logger constructor

	//@}
	
	
	/*!
	 \name Virtual methods : This following methods are overloaded in the derivated classes : XSM_CLOSEST & XSM_MLP & ...
	 */
	//@{
	
	//{
	/// Pure virtual method called to estimates mean cross sections.
	/*!
	 Estimates the mean cross sections evolution according the fresh fuel composition
	 \param IV : Fresh fuel composition
	 \param t : deprecated parameter :
	 \deprecated t : XS update time (used in XSM_Closest)
	 */
	//}
	virtual  EvolutionData GetCrossSections(IsotopicVector IV,double t=0) = 0 ;
	
	/// Check either the IsotopicVector IV is in the validity domain of the models.
	/*!
	 return true if IV is in ValidityDomain
	 return false + a warning if IV is not in ValidityDomain
	 \param IsotopicVector IV, Fresh fuel composition
	 \param t : deprecated parameter :
	 \deprecated t : XS update time (used in XSM_Closest)
	 */
	virtual  bool isIVInDomain(IsotopicVector IV) ;
	//@}
	
	
	
	void SetZAIThreshold(int Z_Threshold){fZAIThreshold = Z_Threshold;}//!< Set the Z threshold : ZAI with Z < fZAIThreshold are not manage by CLASS
	int  GetZAIThreshold(){return fZAIThreshold;}//!< Get the Z threshold

	protected :
	map< ZAI, pair<double,double> > fFreshFuelDomain; //!< Fresh fuel range : map<ZAI<min edge ,max edge >>
	int fZAIThreshold;	//!< Z threshold for handling nuclei mean cross section (take only ZAI reaction of Z>=fZAIThresold)
};

#endif

