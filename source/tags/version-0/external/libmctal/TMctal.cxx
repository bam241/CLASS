#include "TMctal.hxx"
#include <stdio.h>
#include <iomanip>

//________________________________________________________________________
void TMctal::Init() 
{
	fCode[0]=fVersion[0]=fProbId[0]='\0';
	fTallyTab=0;
	fKcode=0;
	fDumps=fNbTal=0;
	fNbPert=0;
	fNPS=fRNR=0;
}

//________________________________________________________________________
unsigned TMctal::Read(istream& s, bool verbose) 
{
	// Read m file
	// Returns 0 if error

	char tmp_str[81];
	unsigned long i;
	unsigned long temp;
	unsigned int champs;

	Free();
	Init();

	if (verbose)
		cout << "\nStart to read m file...\n";
	
	// premiere ligne.

	s.getline(tmp_str,81);
	if (!s.good()) 
	{
		cerr << "Error in TMctal::Read : first line not found\n";
		return 0;
	}

	champs=sscanf(tmp_str,"%8c%8c%19c%ld%lld%lld",fCode,fVersion,fProbId,&fDumps,&fNPS,&fRNR);
	if (champs!=6) 
	{
		cerr << "Error in TMctal::Read : wrong interpretation of the first line\n";
		return 0;
	}

	fCode[8]='\0';
	if (verbose)
		cout << "Code : " << fCode << '\n';
	fVersion[8]='\0';
	if (verbose)
		cout << "Version : " << fVersion << '\n';
	fProbId[19]='\0';
	if (verbose) 
	{
		cout << "Problem ID : " << fProbId << "\n\n";
		cout << fDumps << " dumps, " << fNPS << " source particles.\n";
		cout << fRNR << " random numbers used.\n";
	}

	// comment line

	if (!fCom.Read(s,verbose)) return 0;

	// tally & perturbation number

	s.getline(tmp_str,81);
	if (!s.good()) 
	{
		cerr << "Error in TMctal::Read : tallies & perturbations number not found\n";
		return 0;
	}

	champs=sscanf(tmp_str,"ntal%ldnpert%ld",&fNbTal,&fNbPert);
	if (champs<1){
		cerr << "Error in TMctal::Read : wrong interpretation of tallies & perturbations number\n";
		return 0;
	}
	
	if (verbose)
		cout << fNbTal << " tallies\n";
		
	if (champs==2) 
	{
		if (verbose)
			cout << fNbPert << "perturbations\n";
	}
	else
		fNbPert=1;
	// je ne suis pas sur du rôle de npert.
	// le nombre de tallies à lire est-il ntal ou ntal*npert ?
	// quoi qu'il en soit c'est facile à changer.
	// cf les /**/
	
	//cout << endl;
	
	
	// passer la liste des numéros de tallies

	for (i=0;i<fNbTal*fNbPert;i++) 
	{																				/**/
		s>>temp;
		//s.scan("%ld",&temp);
		if (!s.good()) 
		{
			cerr << "Error in TMctal::Read : read list of tallies failed\n";
			return 0;
		}
	}
	s.getline(tmp_str,81);
	if (!s.good()) 
	{
		cerr << "Error in TMctal::Read : end of tally list\n";
		return 0;
	}


	// creation et lecture de tous les tallies

	fTallyTab = new TMTally[fNbTal*fNbPert];	/**/

	for (i=0;i<fNbTal*fNbPert;i++)				/**/
		if (!fTallyTab[i].Read(s,verbose))
			return 0;

	s.getline(tmp_str,81);
	if (!s.good()) 
	{ // no kcode
		return 1;
	}
	if (tmp_str[0]!='k' || tmp_str[1]!='c' || tmp_str[2]!='o' || tmp_str[3]!='d' || tmp_str[4]!='e') 
	{
		cerr << "Line " << tmp_str << " does not begin by kcode" << endl;
		return 0;
	}

	unsigned long nc,sc,t;
	sscanf(tmp_str+5,"%ld %ld %ld",&nc,&sc,&t);
	fKcode=new TMKcode(nc,sc,t);
	if (!fKcode->Read(s,verbose)) return 0; 
	
	return 1;
}


//________________________________________________________________________
unsigned TMctal::Write(ostream& s) 
{
	unsigned long i;

	if (!fNPS) return 0;

	s << fCode << fVersion << fProbId;
	s<<setw(5)<<fDumps<<setw(11)<<fNPS<<setw(15)<<fRNR<<endl;
	if (!s.good()) return 0;

	if(!fCom.Write(s)) return 0;

	s<<"ntal"<<setw(6)<<fNbTal;
	if(fNbPert>1) /**/
		s<<setw(6)<<fNbPert;
	s << endl;
	if(!s.good()) return 0;

	for (i=0;i<fNbTal;i++)
		 s<<setw(5)<<fTallyTab[i].fNo;
		 s<<endl;
	if(!s.good()) return 0;

	for (i=0;i<fNbTal;i++)
		if(!fTallyTab[i].Write(s)) return 0;

	if (fKcode)
		if (!fKcode->Write(s)) return 0;

	return 1;
}

//________________________________________________________________________
void TMctal::Infos() 
{
	unsigned long i;

	cout << "Code : " << fCode << '\n';
	cout << "Version : " << fVersion << '\n';
	cout << "Problem ID : " << fProbId << "\n\n";
	cout << fDumps << " dumps, " << fNPS << " source particles.\n";
	cout << fRNR << " random numbers used.\n";
	cout << fNbTal << " tally\n";
	if (fNbPert>1)
		cout << fNbPert << "perturbations\n";
	cout << "Tally List:";
	for (i=0;i<fNbTal*fNbPert;i++) 
	{		 /**/
		cout << ' ' << fTallyTab[i].fNo;
	}
	cout << '\n';
}


//________________________________________________________________________
TMTally *TMctal::GetTally(unsigned n) 
{
	// Renvoie le tally d'indice n

	if (fTallyTab && n<fNbTal*fNbPert) /**/
		return fTallyTab+n;

	return 0;

}

//________________________________________________________________________
unsigned TMctal::Add(TMctal *Other) 
{
	// Add Other to this
	// return 0 if error
	// The sum is performed on fNPS, fRNR and tallies
	// but nor on kcodes (is it true?) neither on TFC
	// norme=0 simple sum, 1 if tallies are normalized
	// per source neutron
	
	unsigned i;
	
	if (!Other) return 0;
	if(Other->GetNbTal() != fNbTal*fNbPert) return 0;

	fNPS+=Other->GetNPS();
	fRNR+=Other->GetRNR();
	fDumps+=Other->GetDumps();

	for (i=0;i<fNbTal*fNbPert;i++)	/**/
		if(!fTallyTab[i].Add(Other->GetTally(i)))
			return 0;

	if (fKcode)
		if (!fKcode->Add(Other->GetKcode()))
			return 0;


	return 1;
}

//________________________________________________________________________
TMctal& TMctal::operator *= (const float& factor) 
{
	
	for (unsigned long i=0;i<fNbTal*fNbPert;i++) /**/
		fTallyTab[i]*=factor;

	if (fKcode)
		(*fKcode)*=factor;

	return *this;
}

//________________________________________________________________________
void TMctal::Free() 
{
	delete [] fTallyTab;	
	delete fKcode;
}
