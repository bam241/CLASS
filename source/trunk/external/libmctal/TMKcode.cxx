#include <libmctal/TMKcode.hxx>
#include <cmath>
#include <iomanip>
 

//________________________________________________________________________
TMKcode::TMKcode(unsigned long ActiveCycles,unsigned long InactiveCycles,unsigned long Taille)
{
	fTotalCycles=ActiveCycles;
	fInactiveCycles=InactiveCycles;
	fTaille=Taille;
}

//________________________________________________________________________
TMKcode::~TMKcode() 
{
	for (unsigned long i=0;i<fKcodeCycleTab.size();i++)
		delete fKcodeCycleTab[i];
	
}

//________________________________________________________________________
unsigned TMKcode::Read(istream& s, bool verbose) 
{

	if (verbose) cout << "Reading " << fTotalCycles << " kcode cycles" << endl;

	for (unsigned long i=0;i<fTotalCycles;i++) 
	{
	TMKcodeCycle *c;
		if (fTaille==19)
			c=new TMKcodeCycle19;
	else
			c=new TMKcodeCycle;
		if (!c->Read(s,verbose)) return 0;
		fKcodeCycleTab.push_back(c);

	}
	if (fKcodeCycleTab.size()!=fTotalCycles) 
	{
		cerr << "TMKcode::Read() : Not all cycles could be read" << endl;
		return 0;
	}

	return 1;

}

//________________________________________________________________________
unsigned TMKcode::Write(ostream& s) 
{
	s<<"kcode"<<setw(5)<<fTotalCycles<<setw(5)<<fInactiveCycles<<setw(5)<<fTaille<<endl;
	if (!s.good())
		return 0;

	for (unsigned long i=0;i<fTotalCycles;i++)
		if(!fKcodeCycleTab[i]->Write(s)) return 0;

	return 1;
}

//________________________________________________________________________
unsigned TMKcode::Add(TMKcode *Other) 
{
	if (fInactiveCycles!=Other->fInactiveCycles)
		return 0;
	if (fTaille!=Other->fTaille)
		return 0;
	while (Other->fTotalCycles<fTotalCycles) 
	{
		delete fKcodeCycleTab.back();
		fKcodeCycleTab.pop_back();
		fTotalCycles--;
	}
	for (unsigned long i=0;i<fTotalCycles;i++)
	if (!fKcodeCycleTab[i]->Add(Other->fKcodeCycleTab[i]))
		return 0;
	return 1;
}

//________________________________________________________________________
TMKcode& TMKcode::operator *= (const float& factor) 
{
	for (unsigned long i=0;i<fTotalCycles;i++)
		(*fKcodeCycleTab[i])*=factor;
	return *this;
}


//________________________________________________________________________
unsigned TMKcodeCycle::Read(istream& s, bool verbose) 
{
	s >> fKEffCollision;
	s >> fKEffAbsorption;
	s >> fKEffTL;
	s >> fLifeCollision;
	s >> fLifeAbsorption;
	if (!s) return 0;
	return 1;
}

//________________________________________________________________________
unsigned TMKcodeCycle::Write(ostream& s) 
{

	s<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase);
	s<<" "<<fKEffCollision<<" "<<fKEffAbsorption;
	s<<" "<<fKEffTL<<" "<<fLifeCollision<<" "<<fLifeAbsorption<<endl;
	if (!s) return 0;
	return 1;
}


//________________________________________________________________________
unsigned TMKcodeCycle::Add(TMKcodeCycle *Other) 
{
	fKEffCollision+=Other->fKEffCollision;
	fKEffAbsorption+=Other->fKEffAbsorption;
	fKEffTL+=Other->fKEffTL;
	fLifeCollision+=Other->fLifeCollision;
	fLifeAbsorption+=Other->fLifeAbsorption;
	return 1;
}

//________________________________________________________________________
void TMKcodeCycle::operator *= (const float& factor) 
{
	fKEffCollision*=factor;
	fKEffAbsorption*=factor;
	fKEffTL*=factor;
	fLifeCollision*=factor;
	fLifeAbsorption*=factor;
}

//________________________________________________________________________
unsigned TMKcodeCycle19::Read(istream& s, bool verbose) 
{
	// retrun 0 if error
	s >> fKEffCollision;
	s >> fKEffAbsorption;
	s >> fKEffTL;
	s >> fLifeCollision;
	s >> fLifeAbsorption;
	s >> fKEffCollisionMoy.fVal;
	s >> fKEffCollisionMoy.fErr;
	s >> fKEffAbsorptionMoy.fVal;
	s >> fKEffAbsorptionMoy.fErr;
	s >> fKEffTLMoy.fVal;
	s >> fKEffTLMoy.fErr;
	s >> fKEffAll.fVal; 
	s >> fKEffAll.fErr; 
	s >> fKEffAllSkip.fVal; 
	s >> fKEffAllSkip.fErr; 
	s >> fLifeAll.fVal;
	s >> fLifeAll.fErr; 
	s >> fHistories;
	s >> fKEffAllFM; 

	if (!s)return 0;

	return 1;
}

//________________________________________________________________________
unsigned TMKcodeCycle19::Write(ostream& s) 
{

	s<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase);
	s<<" "<<fKEffCollision<<" "<<fKEffAbsorption;
	s<<" "<<fKEffTL<<" "<<fLifeCollision<<" "<<fLifeAbsorption<<endl;

	s<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase);
	s<<" "<<fKEffCollisionMoy.fVal<<" "<<fKEffCollisionMoy.fErr;
	s<<" "<<fKEffAbsorptionMoy.fVal<<" "<<fKEffAbsorptionMoy.fErr;
	s<<" "<<fKEffTLMoy.fVal<<endl;
	
	s<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase);
	s<<" "<<fKEffTLMoy.fErr<<" "<<fKEffAll.fVal;
	s<<" "<<fKEffAll.fErr<<" "<<fKEffAllSkip.fVal;
	s<<" "<<fKEffAllSkip.fErr<<endl;

	s<<setprecision(5)<<setiosflags(ios::scientific|ios::uppercase);
	s<<" "<<fLifeAll.fVal<<" "<<fLifeAll.fErr;
	s<<" "<<fHistories<<" "<<fKEffAllFM<<endl;
	
	if (!s) return 0;
	return 1;
}

//________________________________________________________________________
unsigned TMKcodeCycle19::Add(TMKcodeCycle *Other) 
{
	//certainly wrong for FOM !
	fKEffCollision+=Other->fKEffCollision;
	fKEffAbsorption+=Other->fKEffAbsorption;
	fKEffTL+=Other->fKEffTL;
	fLifeCollision+=Other->fLifeCollision;
	fLifeAbsorption+=Other->fLifeAbsorption;
	fKEffCollisionMoy+=((TMKcodeCycle19*)Other)->fKEffCollisionMoy;
	fKEffAbsorptionMoy+=((TMKcodeCycle19*)Other)->fKEffAbsorptionMoy;
	fKEffTLMoy+=((TMKcodeCycle19*)Other)->fKEffTLMoy;
	fKEffAll+=((TMKcodeCycle19*)Other)->fKEffAll;
	fKEffAllSkip+=((TMKcodeCycle19*)Other)->fKEffAllSkip;
	fLifeAll+=((TMKcodeCycle19*)Other)->fLifeAll;
	fHistories+=((TMKcodeCycle19*)Other)->fHistories;
	fKEffAllFM+=((TMKcodeCycle19*)Other)->fKEffAllFM;
	return 1;
}

//________________________________________________________________________
void TMKcodeCycle19::operator *= (const float& factor) 
{
	//certainly wrong for FOM !
	fKEffCollision*=factor;
	fKEffAbsorption*=factor;
	fKEffTL*=factor;
	fLifeCollision*=factor;
	fLifeAbsorption*=factor;
	fKEffCollisionMoy*=factor;
	fKEffAbsorptionMoy*=factor;
	fKEffTLMoy*=factor;
	fKEffAll*=factor;
	fKEffAllSkip*=factor;
	fLifeAll*=factor;
	fKEffAllFM*=factor;
}

