#ifndef _CLASSPlotElement_
#define _CLASSPlotElement_

#include "ZAI.hxx"


#include <vector>
#include <iostream>

using namespace std;
	//________________________________________________________________________
	//
	//		CLASSPlotElement
	//@{
	//@}
	//________________________________________________________________________
class CLASSPlotElement
{
	public :
	
	CLASSPlotElement(int treeId, int facilityId, int facylitynumber, int IVNumber, ZAI zai);
	CLASSPlotElement(int treeId, int facilityId, int facylitynumber, int IVNumber, int Z, int A, int I);
	
	virtual ~CLASSPlotElement(){}		//@- destructor
	
	int fTreeId;		// Tree Id
	int fFacilityId;	// FacilityType
				// 0 General
				// 1 reactor
				// 2 Stock
				// 3 Pool
				// 4 FabricationPlant

	int fIVNumber;	// Id of IsotopicVector
			// 0 Incycle
			// 1 CumulativeIN
			// 2 CulumativeOUT

	int fFacylityNumber;	// Id of Facility
				// For General :
				// 0 TOTAL
				// 1 INCYCLE
				// 2 WASTE
				// 3 OUTINCOME
				// 4 REACTOR
				// 5 COOLING
				// 6 STOCK
				// 7 FUELFABRICATION
	ZAI fZAI;		 // ZAI Neeeded
	
};



#endif
