#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "XS/XSM_MLP.hxx"			//Load the include for Neural network cross section predictor
#include "Irradiation/IM_RK4.hxx"		//Load the include for Runge Kutta 4 resolution
#include "Equivalence/EQM_MLP_Kinf.hxx"	//Load the include for Neural Network Equivalence Model (PWRMOX)

TEST ( TestIV, getSize ) {
	const int n=10;

	// génération des données
	IsotopicVector iv;
	for ( unsigned int i = 0 ; i < n ; ++i )
	{
		ZAI z(i,i+1,i+2);
		iv.Add(z,i+3.141592653589);
	}

	EXPECT_EQ( iv.GetZAIQuantity() , n );
}
