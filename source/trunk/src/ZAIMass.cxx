#include "ZAIMass.hxx"

	//________________________________________________________________________
	//
	//		ZAI
	//
	//
	//
	//
	//________________________________________________________________________
	//____________________________InClass Operator____________________________
	//________________________________________________________________________


ZAIMass::ZAIMass()
{
	fZAIMass.insert( pair< ZAI,double >( ZAI(90,232,0), 232000000e-6 ) );
	
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,232,0), 232000000e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,233,0), 233000000e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,234,0), 234000000e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,235,0), 235043929.918e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,236,0), 236000000e-6 ) );	
	fZAIMass.insert( pair< ZAI,double >( ZAI(92,238,0), 238050788.247e-6 ) );

	fZAIMass.insert( pair< ZAI,double >( ZAI(94,238,0), 238049559.894e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,239,0), 239052163.381e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,240,0), 240053813.545e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,241,0), 241056851.456e-6 ) );
	fZAIMass.insert( pair< ZAI,double >( ZAI(94,242,0), 242058742.611e-6 ) );
	
	fZAIMass.insert( pair< ZAI,double >( ZAI(95,241,0), 241056829.144e-6 ) );
	
}


ZAIMass::~ZAIMass()
{
	
}

