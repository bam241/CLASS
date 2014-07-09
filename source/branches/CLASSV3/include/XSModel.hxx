
#ifndef _XSMODEL_HXX
#define _XSMODEL_HXX


/*!
 \file
 \brief Header file for XSMODEL class.
 
 
 @author BLG
 @version 1.0
 */
#include "CLASSHeaders.hxx"


using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a XS Interpolator 
This is the class for method related to XS prediction

see @/Models/XS/include/MLPXSModel.hxx
see @/Models/XS/include/ClosestXSModel.hxx
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

	virtual  EvolutionData GetCrossSections(IsotopicVector IV,double t=0) {return 0;} 



};

#endif

