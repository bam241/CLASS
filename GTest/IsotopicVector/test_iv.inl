#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "XSM_MLP.hxx"			//Load the include for Neural network cross section predictor
#include "IM_RK4.hxx"		    //Load the include for Runge Kutta 4 resolution
#include "EQ_OneParameter.hxx"	//Load the include for Neural Network Equivalence Model (PWRMOX)

TEST ( IV, GetSize ) {
	const int n=10;

	// génération des données
	IsotopicVector iv();
	for ( unsigned int i = 0 ; i < n ; ++i )
	{
		ZAI z(i,i+1,i+2);
		iv.Add(z,i+3.141592653589);
	}

	EXPECT_EQ( iv.GetZAIQuantity() , n );
}
