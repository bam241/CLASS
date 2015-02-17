#ifndef _IMMATRIX_DECAY_
#define _IMMATRIX_DECAY_


/*!
 \file
 \brief Header file for IrradiationModel class.
 
 
 @author BaM
 @version 2.0
 */
#include "IrradiationModel.hxx"
#include "TMatrix.h"


using namespace std;


class CLASSLogger;

//-----------------------------------------------------------------------------//
//! Defines an IrradiationModel based on power series of the exponential of the Bateman matrix

/*!
 Define a IM_Matrix_Decay.
 The aim of these class is to solve numericaly the simplified Bateman equations (taking into account decay only) using the
 development in a power series of the exponential of the Bateman matrix.
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EvolutionData;

class IM_Matrix_Decay : public IrradiationModel
{
	
	public :
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// Default constructor
	
	/*!
	 Make a new IM_Matrix_Decay : */
	IM_Matrix_Decay(IsotopicVector IVList);
	//}
	
	/// Logger constructor
	
	/*!
	 Make a new IM_Matrix_Decay : */
	//param log : Use for the log
	IM_Matrix_Decay(CLASSLogger* log, IsotopicVector IVList);
	//}
	
	//@}
	void 	CleanDecay();

	
	IsotopicVector GetDecay(IsotopicVector Mother_IV, double time);
	
	TMatrixT<double>	ExponentialCalculation(TMatrixT<double> myMatrix);
	void    LoadDecay();				//!< Load the decay properties (HL,BR)

	void	NuclearDataInitialization();		//!< Build Decay matrices & read FpYields if any

	private :
	
	TMatrixT<double>	fExponentialDecayMatrix;	//!< Matrix with half life for each nuclei
	IsotopicVector fIVList;
	
};

#endif

