#include "ZAI.hxx"
#include "CLASSHeaders.hxx"


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

ZAI ZAI::operator=(ZAI IVa)
{
	fZ = IVa.Z();
	fA = IVa.A();
	fI = IVa.I();
	fMass =  IVa.GetMass();
	return *this;
}



ZAI::ZAI()
{
		
	fName="";
	fZ=0;
	fA=0;
	fI=0;
	fMass=0;	
		
}
ZAI::~ZAI()
{
		
}

//________________________________________________________________________
ZAI::ZAI(int Z, int A, int I,bool IsMassSet)
{
		
	fZ=Z;
	fA=A;
	fI=I;

	if(!IsMassSet)
		fMass=0;

	else
	{
		fMass=cZAIMass.GetMass(Z,A);
	}

}

