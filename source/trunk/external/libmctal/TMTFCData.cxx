#include "TMTFCData.hxx"
#include <stdio.h>
#include <iomanip>

//________________________________________________________________________
//
//														TMTFCData
//

//________________________________________________________________________
void TMTFCData::Init() 
{
	fNPS=0;
	fVE=0;
	fM=0;
}

//________________________________________________________________________
unsigned TMTFCData::Read(istream& s, bool verbose) 
{
	char tmp_str[81];

	Init();

	s.getline(tmp_str,81);
	if (!s.good()) 
	{
		cerr << "ERROR in TMTFCData::Read\n";
		return 0;
	}

	sscanf(tmp_str,"%lld%lf%lf",&fNPS,&fVE.fVal,&fVE.fErr);

	if(verbose)
	{
		cout<<fNPS<<"        "<<fVE.fVal<<"        "<<fVE.fErr<<endl;
	}
	fVE.fErr*=fVE.fVal;

	return 1;
}

//________________________________________________________________________
unsigned TMTFCData::Write(ostream& s) 
{
	s<<setw(11)<<fNPS<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fVE.fVal;
	if (fVE.fVal==0)
		s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<0.0;
	else
		s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fVE.fErr/fVE.fVal;
	
	s<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<endl;
	//s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fM<<endl;

	s<<resetiosflags(ios::scientific|ios::uppercase)<<setprecision(6);
	if (!s.good())
		return 0;

	return 1;

}

//________________________________________________________________________
TMTFCData& TMTFCData::operator *= (const float& factor) 
{

	fVE*=factor;
	// What about the FOM ?
	return *this;

}
