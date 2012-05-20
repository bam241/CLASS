#ifndef _TMctal_HXX_
#define _TMctal_HXX_

#include <libmctal/TMTally.hxx>
#include <libmctal/TMComment.hxx>
#include <libmctal/TMKcode.hxx>

#include <iostream>
using namespace std;

//________________________________________________________________________
//@{
// 		A MCNP "m" file.
//
//	This class allows users to interact with "m" file of MCNP : it reads the file and allows user to take its info.
//  Kcode are not completly fully implemented and perturbations are not fully understood.
//	Example: 
// <pre>
//@@		TMctal ficm;
//@@		ifstream MyFile("inpm");
//@@		if(!ficm.Read(MyFile,false))			// read the file with verbose mode off
//@@		{
//@@			cerr<<"\t Error when reading inpm."<<endl;
//@@			exit(0);
//@@		}
//@@		TMTally *TheTally=ficm.GetTally(i);		// get the ith tally of the m file
//@@		cout << "Tally number=" << TheTally->fNo << endl;
//@@		int Nbin=TheTally->fLN;					// get the number of place of the tally (cell or surface bin) 
//@@		for(int i=0; i<Nbin; i++)
//@@		{
//@@			double TallyValue=(TheTally->fVal[i][0][0][0][0][0][0][0]).fVal;
//@@			double TallyRelativeError=(TheTally->fVal[i][0][0][0][0][0][0][0]).fErr
//@@		}
//@@		
//</pre>
//	@see TMTally
//	@author WEC
//	@version 1.0
//@}

class TMctal 
{
 public :
	TMctal() {Init();}							//@- Constructor
	virtual ~TMctal() {Free();} 				//@- Destructor
	unsigned Read(istream& s, bool verbose=1);	//@- Read in "m" file
	unsigned Write(ostream& s); 				//@- write to s
	TMTally *GetTally(unsigned n);				//@- returns the nth tally
	TMKcode *GetKcode() {return fKcode;}		//@- returns the Kcode
	void Infos();								//@- print info
	unsigned long GetNbTal() {return fNbTal*fNbPert;}	//@- returns the number of tally (???)
	unsigned long GetDumps() {return fDumps;}			//@- returns dumps number
	unsigned long long GetNPS() {return fNPS;}			//@- returns particle source number
	void SetNPS(unsigned long long nps) {fNPS=nps;}		//@- set particle source number
	unsigned long long GetRNR() {return fRNR;}			//@- returns number of random number used
	TMComment *GetCom() {return &fCom;}					//@- returns problem description comment

	TMctal& operator *= (const float& factor);
	unsigned Add(TMctal *Other);

 protected :
	void Init();
	void Free();

	char fCode[9];				//@- Code name
	char fVersion[9];			//@- Version
	char fProbId[20];			//@- Date, time, ...
	unsigned long fDumps;		//@- Dumps number
	unsigned long long fNPS;	//@- Particle source number
	unsigned long long fRNR;	//@- Number of random number used
	TMComment fCom;				//@- Problem description comment
	unsigned long fNbTal;		//@- Tally number
	unsigned long fNbPert;		//@- Number of apply perturbation
	TMTally *fTallyTab;			//@- Tally array
	TMKcode *fKcode;			//@- Kcode if exists
								  

};

#endif
