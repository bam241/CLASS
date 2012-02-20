#include "TMTFC.hxx"
#include <iomanip>
#include <stdio.h>



//________________________________________________________________________
void TMTFC::Init() 
{
	fData=0;
	fNb=0;
}

//________________________________________________________________________
unsigned TMTFC::Read(istream& s, bool verbose) 
{
	unsigned long i;

	Free();
	Init();
	
	char tmp_str[81];
	s.getline(tmp_str,81);
	s.getline(tmp_str,81);
 
	char tfc[4];
	sscanf(tmp_str,"%s %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",tfc,&fNb,fIndices,fIndices+1,
		fIndices+2,fIndices+3,fIndices+4,fIndices+5,fIndices+6,fIndices+7);

	if (!s.good()) 
	{
		cerr << "ERROR in TMTFC::Read\n";
		return 0;
	}

	if (verbose) 
	{
		cout << "TFC sur " << fNb << " lignes, indices " << fIndices[0] << ' ' \
			<< fIndices[1] << ' ' << fIndices[2] << ' ' << fIndices[3] << ' ' \
			<< fIndices[4] << ' ' << fIndices[5] << ' ' << fIndices[6] << ' ' \
			<< fIndices[7] << '\n';
	}

	fData=new TMTFCData[fNb];
	
	for (i=0;i<fNb;i++)
		if (!fData[i].Read(s,verbose))
			return 0;

	return 1;
}

//________________________________________________________________________
unsigned TMTFC::Write(ostream& s) 
{
	unsigned long i;

	if (!fNb)
		return 1;

	s<<"tfc"<<setw(5)<<fNb<<setw(8)<<fIndices[0]<<setw(8)<<fIndices[1];
	s<<setw(8)<<fIndices[2]<<setw(8)<<fIndices[3]<<setw(8)<<fIndices[4];
	s<<setw(8)<<fIndices[5]<<setw(8)<<fIndices[6]<<setw(8)<<fIndices[7]<<endl;

	if (!s.good()) return 0;

	for (i=0;i<fNb;i++)
		if(!fData[i].Write(s)) return 0;

	return 1;
}

//________________________________________________________________________
TMTFC& TMTFC::operator *= (const float& factor) 
{
	for (unsigned long i=0;i<fNb;i++)
		fData[i]*=factor;
	
	return *this;
}

unsigned TMTFC::Add(TMTFC *Other) 
{

	// Don't know what to do...
	return 1;

}

//________________________________________________________________________
void TMTFC::Free() 
{
	if(fData)
		delete [] fData;
}
