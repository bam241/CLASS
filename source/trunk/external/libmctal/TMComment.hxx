#ifndef _TMComment_HXX_
#define _TMComment_HXX_

#include <iostream>
using namespace std;

//________________________________________________________________________
//@{
// 		MCNP Comment lines.

//	A MCNP comment line may be on more than one line.
//	If the next line begins by a white spece it is again a comment line.
//
//	@author WEC
//	@version 1.0
//@}

class TMComment 
{

 public :
	TMComment() {Init();}						//@- Constructor
	virtual ~TMComment() {Free();}				//@- Destructor
	unsigned Read(istream& s, bool verbose=1);	//@- Read in "m" file
	unsigned Write(ostream& s);					//@- write the TMComment

	char		fLigne[81];	//@- First comment line
	TMComment	*fNext;		//@- the next ones...

 protected :
	void Init();			//@- Init char[] and pointer
	void Free();			//@- free them

};

#endif
