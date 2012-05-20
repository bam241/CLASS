#ifndef _TMTFC_HXX_
#define _TMTFC_HXX_

#include "libmctal/TMTFCData.hxx"

#include <iostream>
using namespace std;

//________________________________________________________________________
//@{
// 		MCNP Tally fluctuation Chart.
//
//	@author WEC
//	@version 1.0
//@}
class TMTFC 
{

 public:

	TMTFC() {Init();}							//@- Constructor
	virtual ~TMTFC() {Free();}					//@- Destructor
	unsigned Read(istream& s, bool verbose=1);	//@- Read in "m" file
	unsigned Write(ostream& s); 				//@- write to s

	unsigned Add(TMTFC *Other);					//@- not implemeted
	TMTFC& operator *= (const float& factor);

 protected:
	void Init();
	void Free();
	
 public:
	unsigned long fNb;			//@- Entry numbers
	unsigned long fIndices[8];	//@- cell index
	TMTFCData *fData;			//@- the table (TFC)

};

#endif
