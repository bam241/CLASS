#include "TMComment.hxx"


//________________________________________________________________________
void TMComment::Init() 
{
	// init to empty comment
	fNext=0;
	fLigne[0]='\0';
}


//________________________________________________________________________
unsigned TMComment::Read(istream& s, bool verbose) 
{
	// Init from the stream s while lignes begin by space
	// return 0 if error.

	Free();
	Init();

	s.getline(fLigne,81);
	if (!s.good()) 
	{
		cerr << "ERROR in TMComment::Read\n";
		return 0;
	}

	if (verbose) cout << fLigne << '\n';
	if (s.peek()==' ') 
	{
		fNext = new TMComment;
		if (!fNext) 
		{
			cerr << "ERROR in TMComment::Read : memory allocation\n";
			return 0;
		}
		if (!fNext->Read(s,verbose))
			return 0;
	}
	
	return 1;
}


//________________________________________________________________________
unsigned TMComment::Write(ostream& s) 
{
	//write the comment in s (space & comments)
	// returns 0 if error

	if (fLigne[0]) 
	{
		s << fLigne << '\n';
		if (!s.good()) return 0;
	}

	if (fNext)
		return fNext->Write(s);

	return 1;
}


//________________________________________________________________________
void TMComment::Free() 
{
	if(fNext)
		delete fNext;
}
