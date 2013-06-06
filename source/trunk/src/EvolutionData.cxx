#include "EvolutionData.hxx"

/*! \class EvolutionData EvolutionData.hxx "../include/EvolutionData.hxx"
 *
 *  Docs for EvolutionData
 */

#include "IsotopicVector.hxx"
#include "LogFile.hxx"
#include "Defines.hxx"
#include "StringLine.hxx"

#include "TMatrixT.h"
#include <TGraph.h>
#include <TString.h>

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

#pragma link C++ class pair<ZAI,TGraph*>;

using namespace std;
	//________________________________________________________________________
	//
	//		EvolutionData
	//
	//
	//
	//
	//________________________________________________________________________



struct my_tolower
{
	char operator()(char c) const
	{
		return std::tolower(static_cast<unsigned char>(c));
	}
};

	//To Lower Case, convert any string in lower case
string tlc(string data)
{
	transform(data.begin(), data.end(), data.begin(), my_tolower());
	return data;
}

ClassImp(EvolutionData)
	//________________________________________________________________________
EvolutionData operator*(EvolutionData const& evol, double F)
{
	DBGL;
	EvolutionData evoltmp;
	
	map<ZAI ,TGraph* > EvolutionData = evol.GetEvolutionData();
	map<ZAI ,TGraph* >::iterator it;
	for(it = EvolutionData.begin(); it != EvolutionData.end(); it++)
	{
		double X[(*it).second->GetN()];
		double Y[(*it).second->GetN()];
		
		for(int i = 0; i < (*it).second->GetN(); i++)
		{
			double y;
			(*it).second->GetPoint( i, X[i], y );
			Y[i] = y*F;
		}
		evoltmp.NucleiInsert( pair<ZAI, TGraph*> ( (*it).first,new TGraph((*it).second->GetN(), X, Y) ) );
		
	}
	evoltmp.SetPower(evol.GetPower()*F);
	return evoltmp;
	DBGL;
}

	//________________________________________________________________________
EvolutionData operator*(double F, EvolutionData const& evol)
{
	DBGL;
	return evol*F;
	DBGL;
}
	//________________________________________________________________________
EvolutionData operator/(EvolutionData const& evol, double F)
{
	DBGL;
	return evol*(1./F);
	DBGL;
}

	//________________________________________________________________________
	//________________________________________________________________________
EvolutionData::EvolutionData()
{
	DBGL;
	fIsCrossSection = false;
	DBGL;
}

	//________________________________________________________________________
EvolutionData::EvolutionData(LogFile* Log)
{
	DBGL;
	fLog = Log;
	fIsCrossSection = false;
	DBGL;
}

	//________________________________________________________________________
EvolutionData::EvolutionData(LogFile* Log, string DB_file, bool oldread, ZAI zai)
{
	DBGL;
	fLog = Log;
	fIsCrossSection = false;
	fDB_file = DB_file;
	
	if(zai != ZAI(0,0,0))
		AddAsStable(zai);
	else
		ReadDB( fDB_file, oldread);		// Read Evolution Produc DB file name
	
	
	DBGL;
	
}

	//________________________________________________________________________
EvolutionData::~EvolutionData()
{
	DBGL;
	DBGL;
}


bool EvolutionData::NucleiInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	DBGL;
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fEvolutionData.insert( zaitoinsert);
	return IResult.second;
	DBGL;
}

bool EvolutionData::FissionXSInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	DBGL;
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fFissionXS.insert( zaitoinsert);
	return IResult.second;
	DBGL;
}

bool EvolutionData::CaptureXSInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	DBGL;
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fCaptureXS.insert( zaitoinsert);
	return IResult.second;
	DBGL;
}

bool EvolutionData::n2nXSInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	DBGL;
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fn2nXS.insert( zaitoinsert);
	return IResult.second;
	DBGL;
}

	//________________________________________________________________________
void EvolutionData::AddAsStable(ZAI zai)
{
	DBGL;
	double time[2] = {0, (500*365.25*3600*24)};
	double quantity[2] = {1., 1.};
	
	fEvolutionData.insert(pair<ZAI ,TGraph* >(zai, new TGraph(2, time, quantity) ) );
	DBGL;
}

	//________________________________________________________________________
Double_t EvolutionData::Interpolate(double t, TGraph& EvolutionGraph)
{
	DBGL;
	TString fOption;
	return (double)EvolutionGraph.Eval(t,0x0,fOption);
	DBGL;
}

	//________________________________________________________________________
TGraph*	EvolutionData::GetEvolutionTGraph(const ZAI& zai)
{
	DBGL;
	map<ZAI ,TGraph *>::iterator it = GetEvolutionData().find(zai) ;
	
	if ( it != GetEvolutionData().end() )
		return it->second;
	else
		return new TGraph();
	
	DBGL;
}

	//________________________________________________________________________
IsotopicVector	EvolutionData::GetIsotopicVectorAt(double t)
{
	DBGL;
	IsotopicVector IsotopicVectorTmp;
	map<ZAI ,TGraph* >::iterator it;
	for( it = fEvolutionData.begin(); it != fEvolutionData.end(); it++ )
	{
		IsotopicVectorTmp.Add( (*it).first, Interpolate(t, *((*it).second)) );
	}
	
	DBGL;
	return IsotopicVectorTmp;
}

	//________________________________________________________________________
double	EvolutionData::GetGetXSForAt(double t, ZAI zai, int ReactionId)
{
	DBGL;
	map<ZAI ,TGraph* > XSEvol;
	switch(ReactionId)
	{
		case 1: XSEvol = GetFissionXS();
			break;
		case 2: XSEvol = GetCaptureXS();
			break;
		case 3: XSEvol = Getn2nXS();
			break;
		default:cout << "!!Error!! !!!EvolutionData!!! \n Wrong ReactionId !!" << endl;
			fLog->fLog << "!!Error!! !!!EvolutionData!!! \n Wrong ReactionId !!" << endl;
			exit(1);
	}
	
	map<ZAI ,TGraph* >::iterator it = XSEvol.find(zai);
	
	
	if (it == XSEvol.end())
		return 0.;
	else
		return Interpolate(t, *((*it).second) );
	DBGL;
}




	//________________________________________________________________________

EvolutionData EvolutionData::GenerateDBFor(IsotopicVector isotopicvector)
{
	DBGL;
	
	map<ZAI, pair<double, map< ZAI, double > > > ZAIDecay;
	
	{	// TMP
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-3,-3,-3), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// PF
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-2,-2,-2), 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-2,-2,-2), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// 232Th
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(90,232,0), pair<double, map< ZAI, double > > ( 2.37944304000000000e+18 , toAdd ) ) );
	}
	{	// 233U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,233,0), pair<double, map< ZAI, double > > ( 5.02396992000000000e+12, toAdd) ) );
	}
	{	// 234U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,234,0), pair<double, map< ZAI, double > > ( 7.74739080000000000e+12, toAdd) ) );
	}
	{	// 235U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,235,0), pair<double, map< ZAI, double > > ( 2.22165504000000000e+16, toAdd) ) );
	}
	{	// 236U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(90,232,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,236,0), pair<double, map< ZAI, double > > ( 7.39078992000000000e+14, toAdd) ) );
	}
	{	// 238U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,238,0), pair<double, map< ZAI, double > > ( 1.40999356800000000e+17, toAdd) ) );
	}
	{	// 237Np
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(91,233,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(93,237,0), pair<double, map< ZAI, double > > ( 6.76594944000000000e+13, toAdd) ) );
	}
	{	// 238Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,234,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,238,0), pair<double, map< ZAI, double > > ( 2.76760152000000000e+09, toAdd) ) );
	}
	{	// 239Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,235,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,239,0), pair<double, map< ZAI, double > > ( 7.60853736000000000e+11, toAdd) ) );
	}
	{	// 240Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,236,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,240,0), pair<double, map< ZAI, double > > ( 2.07049413600000000e+11, toAdd) ) );
	}
	{	// 241Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,241,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,241,0), pair<double, map< ZAI, double > > ( 4.52062620000000000e+08, toAdd) ) );
	}
	{	// 242Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,238,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,242,0), pair<double, map< ZAI, double > > ( 1.18341000000000000e+13, toAdd) ) );
	}
	{	// 241Am
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,241,0), pair<double, map< ZAI, double > > ( 1.36518177600000000e+10, toAdd) ) );
	}
	{	// 242Am*
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,238,0) , 0.00459) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,0) , 0.99541) );
		
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,242,1), pair<double, map< ZAI, double > > ( 4.44962160000000000e+09, toAdd) ) );
	}
	{	// 243Am
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,243,0), pair<double, map< ZAI, double > > ( 2.32579512000000000e+11, toAdd) ) );
	}
	{	// 242Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,238,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,242,0), pair<double, map< ZAI, double > > ( 1.40659200000000000e+07 , toAdd) ) );
	}
	{	// 243Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 0.9971) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 0.0029) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,243,0), pair<double, map< ZAI, double > > ( 9.18326160000000000e+08, toAdd) ) );
	}
	{	// 244Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,244,0), pair<double, map< ZAI, double > > ( 5.71192560000000000e+08, toAdd) ) );
	}
	{	// 245Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,241,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,245,0), pair<double, map< ZAI, double > > ( 2.65809664800000000e+11, toAdd) ) );
	}
	{	// 246Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,246,0), pair<double, map< ZAI, double > > ( 1.48510065600000000e+11, toAdd) ) );
	}
	{	// 247Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,243,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,247,0), pair<double, map< ZAI, double > > ( 4.92298560000000000e+14, toAdd) ) );
	}
	{	// 248Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,248,0), pair<double, map< ZAI, double > > ( 1.09820448000000000e+13, toAdd) ) );
	}
	
	map<ZAI, map<ZAI, double> > FastDecay;
	{	// 231Th
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(90,231,0), toAdd ) );
	}
	{	// 233Th
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,233,0) , 1) );
		
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(90,233,0), toAdd ) );
	}
	{	// 233Pa
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,233,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(91,233,0), toAdd ) );
	}
	{	// 237U
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(92,237,0), toAdd ) );
	}
	{	// 239U
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(92,239,0), toAdd ) );
	}
	{	// 238Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,238,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,238,0), toAdd ) );
	}
	{	// 239Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,239,0), toAdd ) );
	}
	{	// 240Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,240,0), toAdd ) );
	}
	{	// 241Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,241,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,241,0), toAdd ) );
	}
	{	// 237Pu
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(94,237,0), toAdd ) );
	}
	{	// 243Pu
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(94,243,0), toAdd ) );
	}
	{	// 240Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,240,0), toAdd ) );
	}
	{	// 242Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,242,0) , 0.827) );
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 0.173) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,0), toAdd ) );
	}
	{	// 244Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,244,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,244,0), toAdd ) );
	}
	{	// 245Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,245,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,245,0), toAdd ) );
	}
	{	// 249Cm
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(96,249,0), toAdd ) );
	}
	
	
	map<ZAI, map<ZAI, double> > Capture;
	{	// 241Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,0) , 0.086) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,1) , 0.914) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,241,0), toAdd ) );
	}
	{	// 242Am*
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 1) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,1), toAdd ) );
	}
	map<ZAI, int> index_inver;
	map<int, ZAI> index;
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			index.insert( pair<int, ZAI > ( i, (*it).first ) );
			index_inver.insert( pair<ZAI, int > ( (*it).first , i ));
			i++;
		}
	}
	
	TMatrixT<double> DecayMatrix = TMatrixT<double>(index.size(),index.size());
	for(int i = 0; i < (int)index.size(); i++)
		for(int j = 0; j < (int)index.size(); j++)
			DecayMatrix[i][j] = 0;
	
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			map< ZAI, double >::iterator it2;
			map< ZAI, double > decaylist = (*it).second.second;
			for(it2 = decaylist.begin(); it2!= decaylist.end(); it2++)
			{
				
				map<ZAI, int >::iterator it3 = index_inver.find( (*it2).first );
				if( it3 != index_inver.end() )
					DecayMatrix[(*it3).second][i] = log(2.)/(*it).second.first * (*it2).second;
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = FastDecay.find( (*it2).first );
					
					if( it4 == FastDecay.end() )
					{
						cout << "Problem in FastDecay for nuclei " << (*it2).first.Z() << " " << (*it2).first.A() << " " << (*it2).first.I() << endl;
						exit(1);
					}
					
					map< ZAI, double >::iterator it5;
					map< ZAI, double > decaylist2 = (*it4).second;
					for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
					{
						it3 = index_inver.find( (*it5).first );
						if( it3 == index_inver.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it2).first.Z() << " " << (*it2).first.A() << " " << (*it2).first.I() << endl;
							exit(1);
						}
						DecayMatrix[(*it3).second][i] = log(2.)/(*it).second.first * (*it2).second * (*it5).second;
					}
					
				}
			}
			DecayMatrix[i][i] += -log(2.)/(*it).second.first;
			i++;
			
			
		}
	}
	
	
	vector< TMatrixT<double> > NMatrix ;//  TMatrixT<double>(decayindex.size(),1))
	double NormFactor = 1;
	{
		IsotopicVector WantedHMIV = 	  isotopicvector.GetSpeciesComposition(90)
		+ isotopicvector.GetSpeciesComposition(92)
		+ isotopicvector.GetSpeciesComposition(93)
		+ isotopicvector.GetSpeciesComposition(94)
		+ isotopicvector.GetSpeciesComposition(95)
		+ isotopicvector.GetSpeciesComposition(96);
		
		IsotopicVector DBHMIV =   GetIsotopicVectorAt(0).GetSpeciesComposition(90)
		+ GetIsotopicVectorAt(0).GetSpeciesComposition(92)
		+ GetIsotopicVectorAt(0).GetSpeciesComposition(93)
		+ GetIsotopicVectorAt(0).GetSpeciesComposition(94)
		+ GetIsotopicVectorAt(0).GetSpeciesComposition(95)
		+ GetIsotopicVectorAt(0).GetSpeciesComposition(96);
		
		NormFactor = Norme(WantedHMIV)/ Norme(DBHMIV);
	}
	
	{	// Filling the t=0 State;
		map<ZAI, double > isotopicquantity = isotopicvector.GetIsotopicQuantity();
		TMatrixT<double>  N_0Matrix =  TMatrixT<double>( index.size(),1) ;
		
		map<ZAI, double >::iterator it ;
		for(int i = 0; i < (int)index.size(); i++)
			N_0Matrix[i] = 0;
		
		for(it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{
			
			map<ZAI, int >::iterator it2;
			
			if( (*it).first.Z() < 90 )
				it2 = index_inver.find( ZAI(-2,-2,-2) );
			else it2 = index_inver.find( (*it).first );
			
			if(it2 == index_inver.end() )				//If not in index should be TMP, can't be fast decay for new Fuel !!!
				it2 = index_inver.find( ZAI(-3,-3,-3) );
			
			N_0Matrix[ (*it2).second ][0] = (*it).second ;
			
			
		}
		NMatrix.push_back(N_0Matrix);
	}
	
		//-------------------------//
		//--- Perform Evolution ---//
		//-------------------------//
	double timevector[fEvolutionData.begin()->second->GetN()];
	timevector[0] = 0.;
	
	for(int i = 0; i < fEvolutionData.begin()->second->GetN()-1; i++)
	{
		TMatrixT<double> BatemanMatrix = TMatrixT<double>(index.size(),index.size());
		BatemanMatrix = DecayMatrix ;
		double Flux;
		{
			double x,y;
			fFlux->GetPoint(i, x,y);
			Flux = y;
		}
		map<ZAI ,TGraph* >::iterator it;
			// ----------------  A(n,.) X+Y
		for(it = fFissionXS.begin() ; it != fFissionXS.end(); it++)
		{
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double x,y;
				(*it).second->GetPoint(i, x,y);
				BatemanMatrix[ index_inver.find( (*it).first )->second ][index_inver.find( (*it).first )->second] += -y* 1e-24 *Flux;
				BatemanMatrix[1][ index_inver.find( (*it).first )->second] += 2*y* 1e-24 *Flux;
			}
		}
		
			// ----------------  A(n,.)A+1
		for(it = fCaptureXS.begin() ; it != fCaptureXS.end(); it++)
		{
			
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double x,y;
				(*it).second->GetPoint(i, x, y);
				
				BatemanMatrix[index_inver.find( (*it).first )->second][ index_inver.find( (*it).first )->second ] += -y* 1e-24 *Flux;
				
				map<ZAI, map<ZAI, double> >::iterator it3 = Capture.find( (*it).first );
				
				
				if( it3 == Capture.end() )
				{
					map<ZAI, int >::iterator it6 = index_inver.find( ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
					
					
					if( it6 != index_inver.end() )
					{
						BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux ;
						
					}
					else
					{
						map<ZAI, map<ZAI, double> >::iterator it4 = FastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
						
						if( it4 == FastDecay.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A()+1 << " " << (*it).first.I() << endl;
							exit(1);
						}
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							
							it6 = index_inver.find( (*it5).first );
							if( it6 == index_inver.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A() << " " << (*it).first.I() << endl;
								exit(1);
							}
							
							BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux * (*it5).second;
						}
					}
				}
				else
				{
						//					if( (*it3).first.Z() == 90 && (*it3).first.A() == 232) cout << y* 1e-24 *Flux << endl;
					map<ZAI, double>::iterator it4;
					map<ZAI, double> CaptureList = (*it3).second;
					for(it4 = CaptureList.begin(); it4 != CaptureList.end() ; it4++)
					{
						
						map<ZAI, int >::iterator it6 = index_inver.find( (*it4).first );
						
						
						if( it6 != index_inver.end() )
							BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux * (*it4).second ;
						else
						{
							map<ZAI, map<ZAI, double> >::iterator it7 = FastDecay.find( (*it4).first );
							
							if( it7 == FastDecay.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}
							
							map< ZAI, double >::iterator it5;
							map< ZAI, double > decaylist2 = (*it7).second;
							for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
							{
								
								it6 = index_inver.find( (*it5).first );
								if( it6 == index_inver.end() )
								{
									cout << "Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
									exit(1);
								}
									//if( (*it6).first.Z() == 92 && (*it6).first.A() == 233) cout << y* 1e-24 *Flux * (*it5).second << endl;
								BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y * 1e-24 * Flux * (*it5).second * (*it4).second;
							}
						}
						
					}
				}
				
				
			}
		}
			// ----------------  A(n,2n)A-1
		for(it = fn2nXS.begin() ; it != fn2nXS.end(); it++)
		{
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double x,y;
				(*it).second->GetPoint(i, x,y);
				BatemanMatrix[ index_inver.find( (*it).first )->second ][index_inver.find( (*it).first )->second] += -y* 1e-24 *Flux;
				
				
				map<ZAI, int>::iterator it3 = index_inver.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );
				
				if( it3 != index_inver.end() )
					BatemanMatrix[(*it3).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux;
				else
				{
					
					map<ZAI, map<ZAI, double> >::iterator it4 = FastDecay.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );
					
					if( it4 == FastDecay.end() )
					{
						it3 = index_inver.find( ZAI( -3, -3, -3 ) );
						BatemanMatrix[(*it3).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux;
					}
					else
					{
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							
							it3 = index_inver.find( (*it5).first );
							if( it3 == index_inver.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it4).first.Z() << " " << (*it4).first.A() << " " << (*it4).first.I() << endl;
								exit(1);
							}
							BatemanMatrix[(*it3).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux * (*it5).second ;
						}
					}
				}
			}
		}
		
			// ----------------   Evolution
		TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(index.size(),1);
		
		double TStepMax;
		{
			double x,y;
			fFissionXS.begin()->second->GetPoint(i+1, x,y);
			TStepMax = x - timevector[i];
			timevector[i+1] = x;
		}
		
		BatemanMatrix *= TStepMax;
		
		TMatrixT<double> IdMatrix = TMatrixT<double>(index.size(),index.size());
		for(int j = 0; j < (int)index.size(); j++)
			for(int k = 0; k < (int)index.size(); k++)
			{
				if(k == j)	IdMatrix[j][k] = 1;
				else 		IdMatrix[j][k] = 0;
			}
		
		
		TMatrixT<double> BatemanMatrixDL = TMatrixT<double>(index.size(),index.size());   // Order 0 Term from the DL : Id
		TMatrixT<double> BatemanMatrixDLTermN = TMatrixT<double>(index.size(),index.size());  // Addind it;
		
		
		{
			BatemanMatrixDLTermN = IdMatrix;
			BatemanMatrixDL = BatemanMatrixDLTermN;
			
			
			int j = 1;
			double NormN = 0;
			
			do
			{
				NormN = 0;
				
				TMatrixT<double> BatemanMatrixDLTermtmp = TMatrixT<double>(index.size(),index.size());  // Adding it;
				BatemanMatrixDLTermtmp = BatemanMatrixDLTermN;
				BatemanMatrixDLTermN.Mult(BatemanMatrixDLTermtmp, BatemanMatrix );
				
				BatemanMatrixDLTermN *= 1./j;
				BatemanMatrixDL += BatemanMatrixDLTermN;
				
				NormN = 0;
				for(int m = 0; m < (int)index.size(); m++)
					for(int n = 0; n < (int)index.size(); n++)
						NormN += BatemanMatrixDLTermN[m][n]*BatemanMatrixDLTermN[m][n];
				j++;
			} while ( NormN != 0 || j==2);
			
		}
		
		
		NEvolutionMatrix = BatemanMatrixDL * NMatrix.back() ;
		NMatrix.push_back(NEvolutionMatrix);
	}
	
	
	EvolutionData GeneratedDB = EvolutionData(fLog);
	
	for(int i = 0; i < (int)index.size(); i++)
	{
		double ZAIQuantity[NMatrix.size()];
		for(int j = 0; j < (int)NMatrix.size(); j++)
			ZAIQuantity[j] = (NMatrix[j])[i][0];
		
		GeneratedDB.NucleiInsert(pair<ZAI, TGraph*> (index.find(i)->second, new TGraph(NMatrix.size(), timevector, ZAIQuantity) ) );
	}
	GeneratedDB.SetPower(fPower * NormFactor );
	GeneratedDB.SetFuelType(fFuelType );
	GeneratedDB.SetReactorType(fReactorType );
	GeneratedDB.SetHMMass(fHMMass*NormFactor );
	
	return GeneratedDB;
	DBGL;
}







	//________________________________________________________________________

	//________________________________________________________________________




void EvolutionData::ReadDB(string DBfile, bool oldread)
{
	DBGL;
	if(oldread)
	{
		OldReadDB(DBfile);
		return;
	}

	ifstream DecayDB(DBfile.c_str());	// Open the File
	if(!DecayDB)				//check if file is correctly open
	{
		cout << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
		fLog->fLog << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
	}
	vector<double> vTime;
	
	string line;
	int start = 0;
	
	getline(DecayDB, line);
	if( tlc(StringLine::NextWord(line, start, ' ')) != "time")
	{
		cout << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
		fLog->fLog << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
		vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
	
	fFinalTime = vTime.back();
	double Time[vTime.size()];
	for(int i=0; i < (int)vTime.size();i++)
		Time[i] = vTime[i];
	const int NTimeStep = sizeof(Time)/sizeof(double);

	
	enum { Keff=1, Flux, Inv, XSFis, XSCap, XSn2n };
	
	map<string, int> keyword_map;
	keyword_map["keff"] = Keff;
	keyword_map["flux"] = Flux;
	keyword_map["xsfis"] = XSFis;
	keyword_map["xscap"] = XSCap;
	keyword_map["xsn2n"] = XSn2n;
	keyword_map["inv"] = Inv;

	getline(DecayDB, line);

	do
	{
		start = 0;
		switch (keyword_map[tlc(StringLine::NextWord(line, start, ' '))])
		{
			case Keff:
				ReadKeff(line, Time, NTimeStep);
				break;
				
			case Flux:
				ReadFlux(line, Time, NTimeStep);
				break;

			case Inv:
				ReadInv(line, Time, NTimeStep);
				break;
				
			case XSFis:
				ReadXSFis(line, Time, NTimeStep);
				break;
				
			case XSCap:
				ReadXSCap(line, Time, NTimeStep);
				break;
				
			case XSn2n:
				ReadXSn2n(line, Time, NTimeStep);
				break;
				
				
			default:
				break;
		}
		
		getline(DecayDB, line);

		
	}while ( !DecayDB.eof() )
		
		
		DBGL;
}

void EvolutionData::ReadKeff(string line, double* time, int NTimeStep)
{
	DBGL;

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "keff" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"keff\" not found !" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"keff\" not found !" << endl;
		exit(1);
	}
	
	
	double Keff[NTimeStep];

	int i = 0;
	while(start < (int)line.size() && i < NTimeStep )		// Read the Data
	{
		Keff[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
		i++;
	}
	
	fFlux = new TGraph(NTimeStep, time, Keff);			// Add the TGraph
	
	
	DBGL;
	
}

void EvolutionData::ReadFlux(string line, double* time, int NTimeStep)
{
	DBGL;
	
	int start = 0;
	
	if( tlc(StringLine::NextWord(line, start, ' ')) != "flux" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"flux\" not found !" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"flux\" not found !" << endl;
		exit(1);
	}
	
	double Flux[NTimeStep];
	
	int i = 0;
	while(start < (int)line.size() && i < NTimeStep )		// Read the Data
	{
		Flux[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
		i++;
	}
	
	
	fFlux = new TGraph(NTimeStep, time, Flux);			// Add the TGraph
	
	ReadInfo();							// Read the .info associeted

	DBGL;
}


void	EvolutionData::ReadInv(string line, double* time, int NTimeStep)
{
	DBGL;
	
	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "inv" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"inv\" not found !" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"inv\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double Inv[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			Inv[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}
			// Add the TGraph
		fEvolutionData.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, Inv) ) );
	}
	
	DBGL;
}


void	EvolutionData::ReadXSFis(string line, double* time, int NTimeStep)
{
	DBGL;

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xsfis" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsfis\" not found !" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsfis\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double XSFis[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			XSFis[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}

			// Add the TGraph
		fFissionXS.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, XSFis) ) );
	}
	
	DBGL;
}

void	EvolutionData::ReadXSCap(string line, double* time, int NTimeStep)
{
	DBGL;

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xscap" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xscap\" not found !" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xscap\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double XSCap[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			XSCap[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}

			// Add the TGraph
		fCaptureXS.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, XSCap) ) );
	}
	
	DBGL;
}

void	EvolutionData::ReadXSn2n(string line, double* time, int NTimeStep)
{
	DBGL;

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xsn2n" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsn2n\" not found !" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsn2n\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double XSn2n[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			XSn2n[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}
			// Add the TGraph
		fn2nXS.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, XSn2n) ) );
	}
	
	DBGL;
}



void EvolutionData::ReadInfo()
{
	string InfoDBFile  = fDB_file.erase(fDB_file.size()-3,fDB_file.size());
	InfoDBFile += "Info";
	ifstream InfoDB(InfoDBFile.c_str());				// Open the File
	if(!InfoDB)
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
		fLog->fLog << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
		exit(1);
	}
	
	int start = 0;
	string line;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "reactor")
		fReactorType =  StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if (tlc(StringLine::NextWord(line, start, ' ')) == "fueltype")
		fFuelType =  StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "cycletime")
		fCycleTime =  atof(StringLine::NextWord(line, start, ' ').c_str());;
	getline(InfoDB, line); // Assembly HM Mass
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "constantpower")
		fPower =  atof(StringLine::NextWord(line, start, ' ').c_str());

	
	getline(InfoDB, line); // cutoff
	getline(InfoDB, line); // NUmber of Nuclei
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "normalizationfactor")
	{
		double NormFactor = atof(StringLine::NextWord(line, start, ' ').c_str());
		fPower = fPower * NormFactor;
	}
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "finalheavymetalmass")
		fHMMass =  atof(StringLine::NextWord(line, start, ' ').c_str());
}

void EvolutionData::OldReadDB(string DBfile)
{
	DBGL;
	{
		ifstream DecayDB(DBfile.c_str());							// Open the File
		if(!DecayDB)
		{
			cout << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
			fLog->fLog << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
		}
		vector<double> vTime;
		vector<double> vTimeErr;
		
		string line;
		int start = 0;
		
		getline(DecayDB, line);
		if( StringLine::NextWord(line, start, ' ') != "time")
		{
			cout << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
			fLog->fLog << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
			exit (1);
		}
		
		while(start < (int)line.size())
			vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
		
		fFinalTime = vTime.back();
		double Time[vTime.size()];
		for(int i=0; i < (int)vTime.size();i++)
			Time[i] = vTime[i];
		vector<double> vFlux;
		start = 0;
		getline(DecayDB, line);
		string tmp = StringLine::NextWord(line, start, ' ');
		if ( tmp == "keff"  || tmp == "Keff" )
		{
			vector<double> vKeff;
			while(start < (int)line.size())
				vKeff.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
			
			double Keff[vKeff.size()];
			for(int i=0; i < (int)vKeff.size();i++)
				Keff[i] = vKeff[i];
			
			fKeff = new TGraph(vTime.size(), Time, Keff);
			
			start = 0;
			getline(DecayDB, line);
			if (StringLine::NextWord(line, start, ' ') == "flux")
			{
				
				
				while(start < (int)line.size())
					vFlux.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
				
				double Flux[vFlux.size()];
				for(int i=0; i < (int)vFlux.size();i++)
					Flux[i] = vFlux[i];
				
				fFlux = new TGraph(vTime.size(), Time, Flux);
				
			}
		}
		
		
		do
		{
			
			
			start = 0;
			int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
			
			if(A!=0 && Z!=0)
			{
				double DPQuantity[vTime.size()];
				for(int k = 0; k < (int)vTime.size(); k++ )
					DPQuantity[k] = 0;
				
				
				ZAI zaitmp(Z, A, I);
				int i=0;
				while(start < (int)line.size())
				{
					double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
					DPQuantity[i] = (double)DPQuantityTmp;
					i++;
					
				}
				fEvolutionData.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph((int)vTime.size()-1, Time, DPQuantity) ) );
			}
			
			getline(DecayDB, line);
			if(line == "" || line == "CrossSection" ) break;
		}while (!DecayDB.eof() );
		
		if(line == "CrossSection")
		{
			fIsCrossSection = true;
			getline(DecayDB, line);
			
			if (line == "Fission")
			{
				getline(DecayDB, line);
				
				do
				{
					double DPQuantity[vTime.size()];
					for(int k = 0; k < (int)vTime.size(); k++ )
						DPQuantity[k] = 0;
					
					start = 0;
					int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
					if(A!=0 && Z!=0)
					{
						
						
						ZAI zaitmp(Z, A, I);
						int i=0;
						while(start < (int)line.size())
						{
							long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
							DPQuantity[i] = (double)DPQuantityTmp;
							i++;
							
						}
						fFissionXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
					}
					getline(DecayDB, line);
					if(line == "" || line == "Capture" ) break;
				}while (  !DecayDB.eof() );
			}
			
			if (line == "Capture")
			{
				getline(DecayDB, line); // Nuclei is given with "A Z"
				
				do
				{
					double DPQuantity[vTime.size()];
					for(int k = 0; k < (int)vTime.size(); k++ )
						DPQuantity[k] = 0;
					
					
					start = 0;
					int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
					
					if(A!=0 && Z!=0)
					{
						
						
						ZAI zaitmp(Z, A, I);
						int i=0;
						while(start < (int)line.size())
						{
							long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
							DPQuantity[i] = (double)DPQuantityTmp;
							i++;
							
						}
						fCaptureXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
					}
					getline(DecayDB, line); // Nuclei is given with "A Z"
					if(line == "" || line == "n2n" ) break;
				}while ( !DecayDB.eof() );
				
			}
			
			if (line == "n2n")
			{
				
				getline(DecayDB, line); // Nuclei is given with "A Z"
				
				do
				{
					double DPQuantity[vTime.size()];
					for(int k = 0; k < (int)vTime.size(); k++ )
						DPQuantity[k] = 0;
					
					start = 0;
					int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
					
					if(A!=0 && Z!=0)
					{
						
						
						ZAI zaitmp(Z, A, I);
						int i=0;
						while(start < (int)line.size())
						{
							long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
							DPQuantity[i] = (double)DPQuantityTmp;
							i++;
							
						}
						fn2nXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
					}
					getline(DecayDB, line); // Nuclei is given with "A Z"
					if(line == "" ) break;
					
				}while ( !DecayDB.eof() );
			}
			
		}
		DecayDB.close();
		start = 0;
		
		string InfoDBFile  = DBfile.erase(DBfile.size()-3,DBfile.size());
		InfoDBFile += "info";
		ifstream InfoDB(InfoDBFile.c_str());							// Open the File
		if(!InfoDB)
		{
			fLog->fLog << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
			return;
		}
		
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "Reactor")
			fReactorType =  StringLine::NextWord(line, start, ' ');
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "Fueltype")
			fFuelType =  StringLine::NextWord(line, start, ' ');
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "CycleTime")
			fCycleTime =  atof(StringLine::NextWord(line, start, ' ').c_str());;
		getline(InfoDB, line); // Assembly HM Mass
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "ConstantPower")
			fPower =  atof(StringLine::NextWord(line, start, ' ').c_str());
		getline(InfoDB, line); // cutoff
		getline(InfoDB, line); // NUmber of Nuclei
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "NormalizationFactor")
		{
			double NormFactor = atof(StringLine::NextWord(line, start, ' ').c_str());
			fPower = fPower * NormFactor;
			double Flux[vFlux.size()];
			for(int i=0; i < (int)vFlux.size();i++)
				Flux[i] = vFlux[i];
			
			fFlux = new TGraph(vTime.size()-1, Time, Flux);
		}
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "FinalHeavyMetalMass")
			fHMMass =  atof(StringLine::NextWord(line, start, ' ').c_str());
		
	}
	DBGL;
}





