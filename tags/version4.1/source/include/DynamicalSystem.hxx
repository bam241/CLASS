#ifndef _DynamicalSystem_
#define _DynamicalSystem_ 
/*!
 \file 
  \brief Header file for DynamicalSystem class.
*/ 

#include <math.h>
#include <vector>


using namespace std; 

//-----------------------------------------------------------------------------//


      
//!	DynamicalSystem class solves system of differential equations.  
/*!
//	A DynamicalSystem is a base class that solves system differential equations of 1st order.
//	@f[ \frac{d\vec{Y}}{dt} = A\vec{Y}@f]
//	The differential equations are built in DynamicalSystem::BuildEqns ; the method MUST be
//  defined in the derived classes.
//	In this first version only Runge-Kutta method is implemented, but the aim of this class
//	is to provide also other methods such as CRAM (Chebyshev rational approximation method).
//	
// 
// @author PTO.
// @version 1.0
*/
//________________________________________________________________________

class DynamicalSystem
{
 public :

	DynamicalSystem();								  //!< Normal Constructor 
	DynamicalSystem(const DynamicalSystem & DS);	//!< Copy Constructor 
	virtual ~DynamicalSystem();					   //!< Destructor
	
	/*!
	\name Mains attributes of the DynamicalSystem
	*/
	//@{
	int GetNumberOfEquationSize(){return fNVar;}	//!< return the number of equations.
	void SetNumberOfEquationSize(int n){fNVar = n;}	//!< set the number of equations.
	//@}

	/*!
	\name Runge-Kutta related methods
	Algorithms are taken from Numerical Receipes.
	*/
	//@{

	void SetPrecision(double eps = 1e-5){fPrecision = eps;}	//!< set RK precision to change the integration step
	
	//!	Forbid negative value during integration. 
	/*!
	For some quantities (such as nuclei composition), negative values are forbidden.
	But, due to integration step and very fast variation of the integrated variables
	Runge-Kutta wil produce very small negative value. This method is used to force
	negative value to be zero.
	*/	
	void SetForbidNegativeValue(){fIsNegativeValueAllowed = false;}
	
	//!	Runge Kutta calling method. 
	/*!
	// \param YStart: input : the initial condition Y(t1) ; output the final value Y(t2)
	// \param t1: initial time 
	// \param t2: final time
	*/
	void RungeKutta(double *YStart, double t1, double t2, int EquationNumber);
	
	//!	Builds the equations for integration. 
	/*!
	This method is an abstract method ; it MUST be overwritten by derived classes.
	// \param t: time at which the equations are built 
	// \param Y: array of variable at time t 
	// \param dYdt: ode's variable.
	*/	
	virtual void BuildEqns(double t, double *Y, double *dYdt){}	
		
	//@}
	
	/*!
	\name Miscellaneous methods
	*/
	//@{

	//@}
	
 protected :
	//!	Runge Kutta main method. 
	/*!
	// Call by RungeKutta
	// \param y: initial values to integrate
	// \param dydx: ode's equations (variable is x)
	// \param x: variable of integration
	// \param h: step size for integration
	// \param yout: result after integration
	*/						
	void RK4(double *y, double *dydx, double x, double h, double *yout);	
	//!	Adaptative Step Size method for RK. 
	/*!
	// Call by RK4
	// \param y: initial values to integrate
	// \param dydx: ode's equations (variable is x)
	// \param x: new value of the variable after the adaptative step
	// \param htry: try step size for integration
	// \param eps: precision
	// \param yscal: result after hdid step integration
	// \param hdid: did step size for integration
	// \param hnext: next step size for integration
	*/						
	void AdaptStepSize(double *y, double *dydx, double *x, double htry, double eps, double *yscal, double *hdid, double *hnext); 

	int	fNVar;		 //!< The size of the composition vector and /or number of ZAIs involved. 	 
	double	fPrecision;	//!< Precision of the RungeKutta										 	 
	double	fHestimate;	//!< RK Step estimation. 												 	 
	double	fHmin;		//!< RK minimum Step.													 	 
	double	fMaxHdid;	//!< store the effective RK max step
	double	fMinHdid;	//!< store the effective RK min step
	bool	fIsNegativeValueAllowed; //!< whether or not negative value are physical.
}; 

#endif

