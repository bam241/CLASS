#include "CLASSPlotElement.hxx"


const double AVOGADRO = 6.0221418E+23 ; //2006 CODATA recommended value

using namespace std;
	// for gcc3.2.3 only


CLASSPlotElement::CLASSPlotElement(int treeId, int facilityId, int facylitynumber, int IVNumber, ZAI zai)
{
	
	fTreeId = treeId;
	fFacilityId = facilityId;
	fFacylityNumber = facylitynumber;
	fIVNumber = IVNumber;
	fZAI = zai;
	
}


CLASSPlotElement::CLASSPlotElement(int treeId, int facilityId, int facylitynumber, int IVNumber, int Z, int A, int I)
{
	
	fTreeId = treeId;
	fFacilityId = facilityId;
	fFacylityNumber = facylitynumber;
	fIVNumber = IVNumber;
	fZAI = ZAI(Z,A,I);
	
}

