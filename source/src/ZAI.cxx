#include "ZAI.hxx"
#include "stdlib.h"


//const string DEFAULTDATABASE = "DecayBase.dat";
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
ClassImp(ZAI)

ZAI ZAI::operator = (ZAI IVa)
{
	fZ = IVa.Z();
	fA = IVa.A();
	fI = IVa.I();
	return *this;
}



ZAI::ZAI()
{
		
	fZ = 0;
	fA = 0;
	fI = 0;

}


//________________________________________________________________________
ZAI::ZAI(int Z, int A, int I)
{

	if( Z > A )
	{
		cout << "!!!ERROR!!! " << "[" << __FILE__ << ":" << __FUNCTION__ << "]" << endl;
		cout << "!!!ERROR!!!  Z:" << Z << " is higher than A: " << A << endl;
		cout << "!!!ERROR!!!  CLASS did not manage yet anti-mater!!! Update comming soon !!!"  << endl;
		exit(1);
	}

	fZ = Z;
	fA = A;
	fI = I;

}


ZAI::~ZAI()
{
		
}


