#ifndef _IMMATRIX_
#define _IMMATRIX_


/*!
 \file
 \brief Header file for IrradiationModel class.


 @author BaM
 @version 2.0
 */
#include "IrradiationModel.hxx"


using namespace std;


class CLASSLogger;

//-----------------------------------------------------------------------------//
//! Defines an IrradiationModel based on power series of the exponential of the Bateman matrix

/*!
 Define a IM_Matrix.
 The aim of these class is to solve numericaly the Bateman equations using the
 development in a power series of the exponential of the Bateman matrix.

 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EvolutionData;

class IM_Matrix : public IrradiationModel
{

	public :
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// Default constructor
	
	/*!
	 Make a new IM_Matrix : */
	IM_Matrix();
	//}
	
	/// Logger constructor
	
	/*!
	 Make a new IM_Matrix : */
	//param log : Use for the log
	IM_Matrix(CLASSLogger* log);
	//}
	
	//@}


	/// virtual method called to perform the irradiation calculation using a set of cross section.
	/*!
	 Perform the Irradiation Calcultion using the XSSet data
	 \param IsotopicVector IV isotopic vector to irradiate
	 \param EvolutionData XSSet set of corss section to use to perform the evolution calculation
	 */
	virtual EvolutionData GenerateEvolutionData(IsotopicVector IV, EvolutionData XSSet, double Power, double cycletime);
	//}




	private :

	
	
};

#endif

