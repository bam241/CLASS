#ifndef _TMTFCData_HXX_
#define _TMTFCData_HXX_

#include <libValErr/ValErr.hxx>
#include <fstream>
using namespace std;

//________________________________________________________________________
//@{
// 		One Entry in Tally fluctuation Chart.
//
//	@author WEC
//	@version 1.0
//@}
class TMTFCData 
{
	
 public:
	TMTFCData() {Init();}							//@- Constructor
	virtual ~TMTFCData() {} 						//@- Destructor
	unsigned Read(istream& s, bool verbose=1);		//@- Read in "m" file
	unsigned Write(ostream& s); 					//@- write to s

	TMTFCData& operator *= (const float& factor);

 protected:
	void Init();

	unsigned long long fNPS;	//@- source particle number
	ValErr_t fVE;				//@- Tally, Error
	float fM;					//@- Figure Of Merit 
};

#endif
