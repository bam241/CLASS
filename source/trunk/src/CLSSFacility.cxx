#include "CLSSFacility.hxx"
#include "CLASS.hxx"

using namespace std;

	//________________________________________________________________________
	//
	//		CLSSFacility
	//
	//
	//
	//
	//________________________________________________________________________

ClassImp(CLSSFacility)


CLSSFacility::CLSSFacility()
{
	fParc = 0;
	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = 0;
	fDecayDataBase = 0;
}


//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector CLSSFacility::GetDecay(IsotopicVector isotopicvector, cSecond t)
{

	IsotopicVector IV;

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->Evolution(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}

	return IV;
	
}
