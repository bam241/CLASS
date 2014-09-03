#ifndef _IMMATRIX_HXX
#define _IMMATRIX_HXX


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
/*!
 Define a IM_Matrix.
 The aim of these class is perform the calculation of the evolution of a fuel trough irradiation solving numericaly the Bateman using Matrix.


 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EvolutionData;

class IM_Matrix : public IrradiationModel
{

	public :

	IM_Matrix();
	IM_Matrix(CLASSLogger* log);



	/// virtueal method called to perform the irradiation calculation using a set of cross section.
	/*!
	 Perform the Irradiation Calcultion using the XSSet data
	 \param IsotopicVector IV isotopic vector to irradiate
	 \param EvolutionData XSSet set of corss section to use to perform the evolution calculation
	 */
	EvolutionData GenerateEvolutionData(IsotopicVector IV, EvolutionData XSSet, double Power, double cycletime);
	//}




	private :

	
	
};

#endif

