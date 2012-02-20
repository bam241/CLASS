#include "TMTally.hxx"
#include <stdio.h>
#include <iostream>
#include <iomanip>

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

//________________________________________________________________________
void TMTally::Init() 
{
	fNo=fPart=0;
	fLV=0;
	fCV=fEV=fTV=0;
	fVal=0;
	fLN=fDN=fUN=fSN=fMN=fCN=fEN=fTN=1;
	fUT=fST=fMT=fCT=fET=fTT=1;
	fCC=fEC=fTC=0;
}


//________________________________________________________________________
unsigned TMTally::Read(istream& s, bool verbose) 
{
	// Read tallies in s
	// This is the main part of the "m" file
	//
	// On utilise comme taille de cellule dc * dt * d(ln E)
	// Cas d'un bin total : 1
	// Le premier bin en énergie ne veut rien dire (échelle log)
	//
	// Returns 0 if error

	char tmp_str[81];
	char cle;
	unsigned long i;
	int champs;
	unsigned long l,d,u,sg,m,c,e,t;

	Free();
	Init();

	// ligne de définition

	s.getline(tmp_str,81);
	if (!s.good()) 
	{
		cerr << "Error in TMTally::Read : First line not read\n";
		return 0;
	}

	champs=sscanf(tmp_str,"tally%ld%ld",&fNo,&fPart);
	if (champs!=2) 
	{
		cerr << "Error in TMTally::Read : incorrect interpretation of the first line\n";
		return 0;
	}
	
	if (verbose) 
	{
		cout << "Reading Tally " << fNo << " of ";
		switch(fPart) 
		{
			case 1 :
				cout << "neutrons.\n";
				break;
			case 2 :
				cout << "protons.\n";
				break;
			case 3 :
				cout << "neutrons and protons.\n";
				break;
			case 4 :
				cout << "electrons.\n";
				break;
			case 5 :
				cout << "neutrons and electrons.\n";
				break;
			case 6 :
				cout << "protons and electrons.\n";
				break;
			case 7 :
				cout << "neutrons, protons and electrons.\n";
				break;
			default :
				cout << "Unknown particles !!!\n";
		}
	}

	// comments
	
	if (s.peek()==' ' && !fCom.Read(s,verbose))
		return 0;
	
	do 
	{
		s.getline(tmp_str,81);
		if (!s.good()) 
		{
			cerr << "Error in TMTally::Read : Tally dimension not found\n";
			return 0;
		}

		cle=tmp_str[0];
		switch (cle) 
		{	
			case 'f' : // dimension 1 : The places (cell, surface,...)			

					if (sscanf(tmp_str+2,"%ld",&fLN)!=1) 
					{
						cerr << "Error in TMTally::Read : interpretation of place (cell, surface,...) number\n";
						return 0;
					}
					if (fNo%10 == 5)
						break;
					// la liste n'y est que si ce n'est pas un tally détecteur.
					// pour le type 8, je sais pas ...
			
					fLV=new unsigned long[fLN];
			
					for (i=0;i<fLN;i++) 
					{
						s>>fLV[i];
						if (!s.good()) 
						{
							cerr << "Error in TMTally::Read : read of place list\n";
							return 0;
						}
					}
			
					s.getline(tmp_str,81);
					if (!s.good()) 
					{
						cerr << "Error in TMTally::Read : end of place list\n";
						return 0;
					}
			
					if(fLN>1 && verbose) 
					{
						cout << fLN << " Place (l) :";
						for (i=0;i<fLN;i++)
							cout << ' ' << fLV[i];
						cout << '\n';
					}
					break;
			case 'd' : // dimension 2 : total/direct or flagged/unflagged			
					if (sscanf(tmp_str+2,"%ld",&fDN)!=1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of total/direct or flagged/unflagged\n";
						return 0;
					}
	
					if (fDN>1 && verbose)
						cout << fDN << " total/direct or flagged/unflagged (d)\n";
			
					if (fDN==0) fDN=1;
					break;
			case 'u' : // dimension 3 : user bins
					if (tmp_str[1]=='t')
						fUT=1;
					else fUT=0;
					if (sscanf(tmp_str+2,"%ld",&fUN)!=1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of user bin number\n";
						return 0;
					}
					if (fUN==0) 
					{
						fUN=1;
						fUT=1;
					}

					if (fUN>1 && verbose) 
					{
						cout << fUN << " user bins (u)";
						if (fUT)
							cout << " with one total";
						cout << '\n';
					}
					break;
			case 's' : // dimension 4 : segments
					if (tmp_str[1]=='t')
						fST=1;
					else fST=0;
			
					if (sscanf(tmp_str+2,"%ld",&fSN)!=1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of segment number\n";
						return 0;
					}
			
					if (fSN==0) 
					{
						fSN=1;
						fST=1;
					}

					if (fSN>1 && verbose) 
					{
						cout << fSN << " segments (s)";
						if (fST)
							cout << " with one total";
						cout << '\n';
					}
					break;
			case 'm' : // dimension 5 : multiplicators
					if (tmp_str[1]=='t')
						fMT=1;
					else fMT=0;
			
					if (sscanf(tmp_str+2,"%ld",&fMN)!=1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of multiplicator number\n";
						return 0;
					}
			
					if (fMN==0) 
					{
						fMN=1;
						fMT=1;
					}

					if (fMN>1 && verbose) 
					{
						cout << fMN << " multiplicators (m)";
						if (fMT) cout << " with one total";
						cout << '\n';
					}
					break;
			case 'c' : // dimension 6 : cosines
					if (tmp_str[1]=='t') fCT=1;
					else fCT=0;
			
					champs=sscanf(tmp_str+2,"%ld%ld",&fCN,&fCC);
					if (champs<1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of cosine values\n";
						return 0;
					}

					if (fCN==0) 
					{
						fCN=1;
						fCT=1;
						break;
					}

					if (champs==1)
						fCC=0;

					if (verbose) 
					{
						cout << fCN;
						if (fCC)
							cout << " Center values";
						else
							cout << " bins";
	
						cout << " de cosinus";
	
						if (fCT)
							cout << " with one total";
						cout << '\n';
					}
			
					if (fCC) 
					{
						cerr << "Error in TMTally::Read : Center values not supported\n";
						return 0;
					}

					fCV=new float[fCN-fCT]; /**/

					for (i=0;i<fCN-fCT;i++) 
					{
						s>>fCV[i];
						if (!s.good()) 
						{
							cerr << "Error in TMTally::Read : incorrect list of des cosinus\n";
							return 0;
						}
					}
					s.getline(tmp_str,81);
					if (!s.good()) 
					{
						cerr << "Error in TMTally::Read : wrong end list of cosines\n";
						return 0;
					}
					break;
			case 't' : // dimension 6 : times		
					if (tmp_str[1]=='t') fTT=1;
					else fTT=0;
			
					champs=sscanf(tmp_str+2,"%ld%ld",&fTN,&fTC);
					if (champs<1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of time values\n";
						return 0;
					}

					if (fTN==0) 
					{
						fTN=1;
						fTT=1;
						break;
					}

					if (champs==1)
						fTC=0;
			
					if (verbose) 
					{
						cout << fTN;
						if (fTC)
							cout << " center values";
						else
							cout << " bins";
						cout << " de temps";
	
						if (fTT)
							cout << " with one total";
						cout << '\n';
					}
			
					if (fTC) 
					{
						cerr << "Error in TMTally::Read : Center values not supported\n";
						return 0;
					}
					fTV=new float[fTN-fTT];
				
					for (i=0;i<fTN-fTT;i++) 
					{
						s>>fTV[i];
						if (!s.good()) 
						{
							cerr << "Error in TMTally::Read : incorrect list of time\n";
							return 0;
						}
					}

					s.getline(tmp_str,81);
					if (!s.good()) 
					{
						cerr << "Error in TMTally::Read : wrong end list of times\n";
						return 0;
					}
					break;
			case 'e' : // dimension 7 : energies
					if (tmp_str[1]=='t')
			 			fET=1;
					else fET=0;
			
					champs=sscanf(tmp_str+2,"%ld%ld",&fEN,&fEC);
					if (champs<1) 
					{
						cerr << "Error in TMTally::Read : incorrect interpretation of energy number\n";
						return 0;
					}

					if (fEN==0) 
					{
						fEN=1;
						fET=1;
						break;
					}

					if (champs==1)fEC=0;

					if (verbose) 
					{
						cout << fEN;
						if (fEC)
							cout << " center values";
						else
							cout << " bins";
						cout << " of energies";

						if (fET)
							cout << " with one total";
						cout << '\n';
					}
					if (fEC) 
					{
						cerr << "Error in TMTally::Read : Center values not supported\n";
						return 0;
					}
					
					fEV=new float[fEN-fET];
					for (i=0;i<fEN-fET;i++) 
					{
						s>>fEV[i];
						if (!s.good()) 
						{
							cerr << "Error in TMTally::Read : incorrect list of energies\n";
							return 0;
						}
					}
					s.getline(tmp_str,81);
					if (!s.good()) 
					{
						cerr << "Error in TMTally::Read : wrong end list of energies\n";
						return 0;
					}
					break;
			case 'v' : break;
			default :
				cerr << "Error in TMTally::Read : Undefined key : " << cle << '\n';
				return 0;
			
		}
	} while (cle!='v');
	
	fVal=new ValErr_t*******[fLN];
	if(!fVal) 
	{
		cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
		return 0;
	}
	for (l=0;l<fLN;l++) 
	{
		fVal[l]=new ValErr_t******[fDN];
		if(!fVal[l]) 
		{
			cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
			return 0;
		}
		for (d=0;d<fDN;d++) 
		{
			fVal[l][d]=new ValErr_t*****[fUN];
			if(!fVal[l][d]) 
			{
				cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
				return 0;
			}
			for (u=0;u<fUN;u++) 
			{
				fVal[l][d][u]=new ValErr_t****[fSN];
				if(!fVal[l][d][u]) 
				{
					cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
					return 0;
				}
				for (sg=0;sg<fSN;sg++) 
				{
					fVal[l][d][u][sg]=new ValErr_t***[fMN];
					if(!fVal[l][d][u][sg]) 
					{
						cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
						return 0;
					}
					for (m=0;m<fMN;m++) 
					{
						fVal[l][d][u][sg][m]=new ValErr_t**[fCN];
						if(!fVal[l][d][u][sg][m]) 
						{
							cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
							return 0;
						}
						for (c=0;c<fCN;c++) 
						{
							fVal[l][d][u][sg][m][c]=new ValErr_t*[fEN];
							if(!fVal[l][d][u][sg][m][c]) 
							{
								cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
								return 0;
							}
							for (e=0;e<fEN;e++) 
							{
								fVal[l][d][u][sg][m][c][e]=new ValErr_t[fTN];
								if(!fVal[l][d][u][sg][m][c][e]) 
								{
									cerr << "Erreur à l'allocation mémoire pour le contenu du tally\n";
									return 0;
								}
								for (t=0;t<fTN;t++) 
								{
									s>>fVal[l][d][u][sg][m][c][e][t].fVal>>fVal[l][d][u][sg][m][c][e][t].fErr;
									if (!s.good()) 
									{
										cerr << "Error in TMTally::Read : wrong tally contains\n";
										return 0;
									}
									fVal[l][d][u][sg][m][c][e][t].fErr*=fVal[l][d][u][sg][m][c][e][t].fVal;						
								}
							}
						}
					}
				}
			}
		}
	}

	if (!fTFC.Read(s,verbose)) return 0;

	return 1;
}

//________________________________________________________________________
unsigned TMTally::Write(ostream& s) 
{
	// Write a TMTally in s
	// No beautifull appearence but compatible with TMTally::Read
	// Returns 0 if error

	unsigned long i;
	unsigned long l,d,u,sg,m,c,e,t;

	if (!fNo) return 0;

	s<<"tally"<<setw(5)<<fNo<<setw(5)<<fPart<<endl;
	if (!s.good()) return 0;

	if(!fCom.Write(s)) return 0;

	s<<"f"<<setw(9)<<fLN<<endl;
	if(fNo%10 != 5)
		for (i=0;i<fLN;i++)
			s<<setw(7)<<fLV[i];
	s << '\n';
	if (!s.good()) return 0;

	s<<"d"<<setw(9)<<fDN<<endl;
	if (!s.good()) return 0;

	if (fUN==1 && fUT==1)
		s << "u        0\n";
	else 
	{
		s << "u";
		if (fUT) 
		{
			s << "t"<<setw(8)<<fUN<<endl;
		}
		else
			 s<<setw(9)<<fUN<<endl;
	}
	if (!s.good()) return 0;
	
	if (fSN==1 && fST==1)
		s << "s        0\n";
	else 
	{
		s << "s";
		if (fST) 
		{
			s << "t"<<setw(8)<<fSN<<endl;
		}
		else
			s<<setw(9)<<fSN<<endl;
	}
	if (!s.good()) return 0;

	if (fMN==1 && fMT==1)
		s << "m        0\n";
	else 
	{
		s << "m";
		if (fMT) 
		{
			s << "t"<<setw(8)<<fMN<<endl;
		}
		else
			s<<setw(9)<<fMN<<endl;
	}
	if (!s.good()) return 0;

	if (fCN==1 && fCT==1)
		s << "c        0\n";
	else 
	{
		s << "c";
		if (fCT) 
		{
			s << "t"<<setw(8)<<fCN<<endl;
		}
		else
			s<<setw(9)<<fCN<<endl;
		for (c=0;c<fCN-fCT;c++) 
		{
			s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fCV[c];
			if (c%6==5 || c==fCN-fCT-1)
				s << endl;
		}
	}
	if (!s.good()) return 0;

	if (fEN==1 && fET==1)
		s << "e        0\n";
	else 
	{
		s << "e";
		if (fET) 
		{
			s << "t"<<setw(8)<<fEN<<endl;
		}
		else
			s<<setw(9)<<fEN<<endl;
		for (e=0;e<fEN-fET;e++) 
		{
			s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fEV[e];
			if (e%6==5 || e==fEN-fET-1)
				s << endl;
		}
	}
	if (!s.good()) return 0;

	if (fTN==1 && fTT==1)
		s << "t        0\n";
	else 
	{
		s << "t";
		if (fTT) 
		{
			s << "t"<<setw(8)<<fTN<<endl;
		}
		else
			s<<setw(9)<<fTN<<endl;;
		for (t=0;t<fTN-fTT;t++) 
		{
			s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fTV[t];
			if (t%6==5 || t==fTN-fTT-1)
				s << endl;
		}
	}
	if (!s.good()) return 0;

	s << "vals\n";
	unsigned long long indice=0;
	for (l=0;l<fLN;l++)
		for (d=0;d<fDN;d++)
			for (u=0;u<fUN;u++)
				for (sg=0;sg<fSN;sg++)
					for (m=0;m<fMN;m++)
						for (c=0;c<fCN;c++)
							for (e=0;e<fEN;e++)
								for (t=0;t<fTN;t++) 
								{
									s<<"        "<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase)<<fVal[l][d][u][sg][m][c][e][t].fVal<<" ";
									if (fVal[l][d][u][sg][m][c][e][t].fVal!=0)
										s<<resetiosflags(ios::scientific|ios::uppercase)<<setiosflags(ios::fixed)<<setprecision(4)<<fVal[l][d][u][sg][m][c][e][t].fErr/fVal[l][d][u][sg][m][c][e][t].fVal<<resetiosflags(ios::fixed);
									else
										s<<resetiosflags(ios::scientific|ios::uppercase)<<setiosflags(ios::fixed)<<setprecision(4)<<0.0<<resetiosflags(ios::fixed);
									if (!s.good()) return 0;
									if (indice%4==3 || indice==fLN*fDN*fUN*fSN*fMN*fCN*fEN*fTN-1)
										s << endl;
									indice++;
								}

	s<<resetiosflags(ios::scientific|ios::uppercase)<<resetiosflags(ios::fixed)<<setprecision(6);
	if (!fTFC.Write(s)) return 0;

	return 1;
}

//________________________________________________________________________
unsigned TMTally::Add(TMTally *Other) 
{
	// Add  Other Tally to this
	// returns 0 if error
	
	unsigned long l,d,u,sg,m,c,e,t;
	
	if(!Other) return 0;
	if(fNo!=Other->fNo) return 0;
	if(fPart!=Other->fPart) return 0;
	if(fLN!=Other->fLN) return 0;
	if(fNo%10 != 5)
		for(l=0;l<fLN;l++)
			if(fLV[l]!=Other->fLV[l]) return 0;
	if(fDN!=Other->fDN) return 0;
	if(fUN!=Other->fUN) return 0;
	if(fUT!=Other->fUT) return 0;
	if(fSN!=Other->fSN) return 0;
	if(fST!=Other->fST) return 0;
	if(fMN!=Other->fMN) return 0;
	if(fMT!=Other->fMT) return 0;
	if(fCN!=Other->fCN) return 0;
	if(fCT!=Other->fCT) return 0;
	for(c=0;c<fCN-fCT;c++)
		if(fCV[c]!=Other->fCV[c]) return 0;
	if(fCC!=Other->fCC) return 0;
	if(fEN!=Other->fEN) return 0;
	if(fET!=Other->fET) return 0;
	for(e=0;e<fEN-fET;e++)
		if(fEV[e]!=Other->fEV[e]) return 0;
	if(fEC!=Other->fEC) return 0;
	if(fTN!=Other->fTN) return 0;
	if(fTT!=Other->fTT) return 0;
	for(t=0;t<fTN-fTT;t++)
		if(fTV[t]!=Other->fTV[t]) return 0;
	if(fTC!=Other->fTC) return 0;

	for (l=0;l<fLN;l++)
		for (d=0;d<fDN;d++)
			for (u=0;u<fUN;u++)
				for (sg=0;sg<fSN;sg++)
					for (m=0;m<fMN;m++)
						for (c=0;c<fCN;c++)
							for (e=0;e<fEN;e++)
								for (t=0;t<fTN;t++)
									fVal[l][d][u][sg][m][c][e][t]+=Other->fVal[l][d][u][sg][m][c][e][t];
	
	if(!fTFC.Add(&(Other->fTFC))) return 0;

	return 1;
}


//________________________________________________________________________
void TMTally::Infos() 
{
	// Affiche dans cout des infos concernant le tally

	unsigned long i;

	cout << "\nTally " << fNo << " of ";
	switch(fPart) 
	{
		case 1 :
			cout << "neutrons.\n";
			break;
		case 2 :
			cout << "protons.\n";
			break;
		case 3 :
			cout << "neutrons and protons.\n";
			break;
		case 4 :
			cout << "electrons.\n";
			break;
		case 5 :
			cout << "neutrons and electrons.\n";
			break;
		case 6 :
			cout << "protons and electrons.\n";
			break;
		case 7 :
			cout << "neutrons, protons and electrons.\n";
			break;
		default :
			cout << "Unknown particles !!!\n";
	}

	cout << fLN << " Places (l) :";
	for (i=0;i<fLN;i++)
		cout << ' ' << fLV[i];
	cout << '\n';
	
	if (fDN>1)
		cout << fDN << " total/direct or flagged/unflagged (d)\n";

	if (fUN>1) 
	{
		cout << fUN << " user bins (u)";
		if (fUT)
			cout << " with one total";
		cout << '\n';
	}

	if (fSN>1) 
	{
		cout << fSN << " segments (s)";
		if (fST)
			cout << " with one total";
		cout << '\n';
	}

	if (fMN>1) 
	{
		cout << fMN << " multiplicators (m)";
		if (fMT)
			cout << " with one total";
		cout << '\n';
	}

	if (fCN-fCT != 0) 
	{
		cout << fCN;
		if (fCC)
			cout << " center values";
		else
			cout << " bins";
		cout << " of cosines";
		if (fCT)
			cout << " with one total";
		cout << '\n';
	}

	if (fTN-fTT != 0) 
	{
		cout << fTN;
		if (fTC)
			cout << " center values";
		else
			cout << " bins";
		cout << " of times";
		if (fTT)
			cout << " with one total";
		cout << '\n';
	}

	if (fEN-fET != 0) 
	{
		cout << fEN;
		if (fEC)
			cout << " center values";
		else
			cout << " bins";
		cout << " of energies";
		if (fET)
			cout << " with one total";
		cout << "\n\n";
	}
	
}


//________________________________________________________________________
void TMTally::Free() 
{
	unsigned long l,d,u,s,m,c,e;

	delete [] fLV;
	delete [] fCV;
	delete [] fEV;
	delete [] fTV;

	if (fVal) 
	{
		for (l=0;l<fLN;l++) 
		{
			for (d=0;d<fDN;d++) 
			{
				for (u=0;u<fUN;u++) 
				{
					for (s=0;s<fSN;s++) 
					{
						for (m=0;m<fMN;m++) 
						{
							for (c=0;c<fCN;c++) 
							{
								for (e=0;e<fEN;e++)
									delete [] fVal[l][d][u][s][m][c][e];
								delete [] fVal[l][d][u][s][m][c];
							}
							delete [] fVal[l][d][u][s][m];
						}
						delete [] fVal[l][d][u][s];
					}
					delete [] fVal[l][d][u];
				}
				delete [] fVal[l][d];
			}
			delete [] fVal[l];
		}
		delete [] fVal;
	}
}



//________________________________________________________________________
TMTally& TMTally::operator *= (const float& factor) 
{
	for (unsigned long l=0;l<fLN;l++)
		for (unsigned long d=0;d<fDN;d++)
			for (unsigned long u=0;u<fUN;u++)
				for (unsigned long sg=0;sg<fSN;sg++)
					for (unsigned long m=0;m<fMN;m++)
						for (unsigned long c=0;c<fCN;c++)
							for (unsigned long e=0;e<fEN;e++)
								for (unsigned long t=0;t<fTN;t++)
									fVal[l][d][u][sg][m][c][e][t]*=factor;
	fTFC*=factor;
	return *this;
}


/*
unsigned TMTally::ReBinT(float *TV, unsigned long TN,unsigned long TT) 
{
// Rebinne le tally en temps (eau qui tombe verticalement dans de nouvelles cases).
// Ne marche pas pour des valeurs centrees. Retourne 0 si rate.
// Gere l'erreur.
// Suppose les valeurs de t rangees par ordre croissant.

unsigned long l,d,u,s,m,c,e;
unsigned long i,j;
float deb,fin,dt;

TMVal *NewVals;

if (!fVal)
return 0;
if (fTC)
return 0;

for (l=0;l<fLN;l++) 
{
for (d=0;d<fDN;d++) 
{
for (u=0;u<fUN;u++) 
{
for (s=0;s<fSN;s++) 
{
for (m=0;m<fMN;m++) 
{
for (c=0;c<fCN;c++) 
{
for (e=0;e<fEN;e++) 
{

NewVals=new TMVal[TN];
if(!NewVals) 
{
cout << "Erreur à l'allocation mémoire pour le contenu du tally\n";
return 0;
}

i=j=0;
deb=fin=0;
dt=fTV[0];

NewVals[0].SetValeur(0);
NewVals[0].SetErreur(0);

while (i<TN-TT && j<fTN-fTT) 
{
deb=fin;
fin=fTV[j]<TV[i]?fTV[j]:TV[i];

if (fin>deb)
NewVals[i].Add(&fVal[l][d][u][s][m][c][e][j],1,(fin-deb)/dt);

if (fTV[j]<TV[i]) 
{
j++;
dt=fTV[j]-fTV[j-1];
}
else {
i++;
NewVals[i].SetValeur(0);
NewVals[i].SetErreur(0);
}

}

while(++i<TN-TT) 
{
NewVals[i].SetValeur(0);
NewVals[i].SetErreur(0);
}

if (TT) 
{
NewVals[TN-1].SetValeur(0);
NewVals[TN-1].SetErreur(0);
for (i=0;i<TN-1;i++)
NewVals[TN-1].Add(&NewVals[i],1,1);
}

delete[] fVal[l][d][u][s][m][c][e];
fVal[l][d][u][s][m][c][e]=NewVals;

}
}
}
}
}
}
}

delete fTV;
fTV=new float[TN-TT];
if (!fTV)
return 0;
for (i=0;i<TN-TT;i++)
fTV[i]=TV[i];
fTN=TN;
fTT=TT;

return 1;

}
*/
