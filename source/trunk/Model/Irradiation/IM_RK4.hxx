#ifndef _IMRK4_HXX
#define _IMRK4_HXX


/*!
 \file
 \brief Header file for IrradiationModel class.


 @author BaM
 @version 2.0
 */
#include "DynamicalSystem.hxx"
#include "IrradiationModel.hxx"


using namespace std;


class CLASSLogger;

//-----------------------------------------------------------------------------//
/*!
 Define a IM_RK4
 
 The aim of these class is to solve numericaly the Bateman equations using 
 Runge Kutta 4th order method
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EvolutionData;

class IM_RK4 : public IrradiationModel, DynamicalSystem
{

	public :
	
	/*!
	 \name Constructor
	 */
	//@{
	
	//{
	/// Default constructor
	//
	/*!
	 Make a new IM_Matrix : */
	IM_RK4();
	//}
	
	//{
	/// Logger constructor
	//
	/*!
	 Make a new IM_Matrix :
	\param log : Use for the log
	*/
	IM_RK4(CLASSLogger* Log);
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

	//********* RK4 Method *********//

	//@}
	/*!
	 \name RK4 Method
	 */
	//@{


	using	DynamicalSystem::RungeKutta;
	//!	Pre-treatment Runge-Kutta method.
	/*!
	 // This method does initialization and then call DynamicalSystem::RungeKutta
	 // \param t1: initial time
	 // \param t2: final time
	 */


   	void BuildEqns(double t, double *N, double *dNdt);
	void SetTheMatrixToZero();				//!< Initialize the evolution Matrix
	void ResetTheMatrix();
	void SetTheMatrix(TMatrixT<double> BatemanMatrix);	//!< Set the Evolution Matrix (Bateman equations)
	TMatrixT<double> GetTheMatrix();			//!< return the Evolution Matrix (Bateman equations)

	void SetTheNucleiVectorToZero();			//!< Initialize the evolution Matrix
	void ResetTheNucleiVector();
	void SetTheNucleiVector(TMatrixT<double> NEvolutionMatrix);	//!< Set the Evolution Matrix (Bateman equations)
	TMatrixT<double> GetTheNucleiVector();			//!< return the Evolution Matrix (Bateman equations)
	//@}



	private :

	double	*fTheNucleiVector;	//!< The evolving atoms copied from Material proportions.
	double 	**fTheMatrix;  		//!< The evolution Matrix

	double	fPrecision;	//!< Precision of the RungeKutta
	double	fHestimate;	//!< RK Step estimation.
	double	fHmin;		//!< RK minimum Step.
	double	fMaxHdid;	//!< store the effective RK max step
	double	fMinHdid;	//!< store the effective RK min step
	bool	fIsNegativeValueAllowed; //!< whether or not negative value are physical.


};

#endif

