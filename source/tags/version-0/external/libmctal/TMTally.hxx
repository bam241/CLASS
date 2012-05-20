#ifndef _TMTally_HXX_
#define _TMTally_HXX_

#include <libValErr/ValErr.hxx>
#include <libmctal/TMComment.hxx>
#include <libmctal/TMTFC.hxx>

#include <iostream>

//________________________________________________________________________
//@{
// 		A MCNP Tally.
//
//	It contains its fluctuation table (TFC). Units are MCNP units (MeV, shakes,...). 
//
//  For binning array, If values of bin are center, there is no need to divid by bin width to obtain the mean
//	whereas it should be done if bins are given by there upper bound (default).
//
//	@author WEC
//	@version 1.0
//@}

class TMTally 
{

 public:
	TMTally() {Init();}							//@- Constructor
	virtual ~TMTally() {Free();}				//@- Destructor
	unsigned Read(istream& s, bool verbose=1);	//@- Read in "m" file
	unsigned Write(ostream& s); 				//@- write to s
	void Infos();								//@- print info

	unsigned long fNo;		//@- Tally Number (e.g. 14)
	unsigned long fPart;	//@- Particle type (1=N,2=P,3=NP,4=E,5=NE,6=PE,7=NPE)
	
	unsigned long fLN;		//@- Cell or Surface or detector bin number
	unsigned long *fLV;		//@- Array of their numbers (except for detector)

	unsigned long fDN;		//@- Total/direct or flagged/unflagged
	unsigned long fUN;		//@- Number of user bin
	unsigned char fUT;		//@- Is there one of those for the total (y=1, no=0)

	unsigned long fSN;		//@- Number of bin segment
	unsigned char fST;		//@- Total?

	unsigned long fMN;		//@- Number of multiplicators
	unsigned char fMT;		//@- Total?

	unsigned long fCN;		//@- Cosine bin number
	unsigned char fCT;		//@- Total?
	unsigned long fCC;		//@- Are the values center (1) or upper bound (0) 
	float *fCV;				//@- Array of cosines

	unsigned long fEN;		//@- Number of Energy bins
	unsigned char fET;		//@- Total?
	unsigned long fEC;		//@- Are the values center (1) or upper bound (0) 
	float *fEV;				//@- Array of energies

	unsigned long fTN;		//@- Number of Time bins
	unsigned char fTT;		//@- Total?
	unsigned long fTC;		//@- Are the values center (1) or upper bound (0)
	float *fTV;				//@- Array of times
	//@{
	// 	Tally contains.
	//
	//	Tally contains without division by dE,dt, ...
	//
	// fVal[l][d][u][s][m][c][e][t] where l=cell, d=flag, u=user, s=segment, m=multiplicator
	// 								c=cosine, e=energy, t=time
	//@}	
	ValErr_t ********fVal;
	TMTFC fTFC;				//@- Tally Fluctuation Chart
	TMComment fCom; 		//@- Comment

	TMTally& operator *= (const float& factor);	//@- Multiply a Tally by factor
	unsigned Add(TMTally *Other);				//@- Add 2 tallies
	
	// To be done
	//unsigned Add(TMTally *OtherTally,float Coef1,float Coef2);
	//To be suppress
	//	unsigned ReBinT(float *TV, unsigned long TN, unsigned long TT);

 protected:
	void Init();
	void Free();

};

#endif
