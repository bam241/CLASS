#include "IsotopicVector.hxx"

#include "LogFile.hxx"
#include "Defines.hxx"


#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

//________________________________________________________________________
//________________________________________________________________________
//
//
//
//				IsotopicVector
//
//
//________________________________________________________________________
//________________________________________________________________________




//________________________________________________________________________
//__________________________Operator Overlaoding__________________________
//________________________________________________________________________



//____________________________General Operator____________________________
//________________________________________________________________________

ClassImp(IsotopicVector)

double 	Norme(IsotopicVector IV1,int DistanceType, IsotopicVector DistanceParameter)
{
DBGL;
	IsotopicVector IV; 
	
	return Distance(IV1, IV, DistanceType,DistanceParameter);
DBGL;
}
double DistanceStandard(IsotopicVector IV1, IsotopicVector IV2)
{
DBGL;
	double d2=0;
	IsotopicVector IVtmp = IV1 + IV2;
	map<ZAI ,double> IVtmpIsotopicQuantity = IVtmp.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		double Z1 = IV1.GetZAIIsotopicQuantity( (*it).first );
		double Z2 = IV2.GetZAIIsotopicQuantity( (*it).first );
		d2 += pow(Z1-Z2 , 2 );
	}
	return sqrt(d2);
DBGL;
}
double DistanceAdjusted(IsotopicVector IV1, IsotopicVector IV2, IsotopicVector DistanceParameter)
{
DBGL;
	double d2=0;
	IsotopicVector IVtmp = IV1 + IV2;
	map<ZAI ,double> IVtmpIsotopicQuantity = IVtmp.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		double Z1 = IV1.GetZAIIsotopicQuantity( (*it).first );
		double Z2 = IV2.GetZAIIsotopicQuantity( (*it).first );
		double lambda = DistanceParameter.GetZAIIsotopicQuantity( (*it).first );
		d2 += lambda*abs(Z1-Z2);
	}
	return d2;
DBGL;
}



double Distance(IsotopicVector IV1, IsotopicVector IV2 ,int DistanceType, IsotopicVector DistanceParameter)
{
DBGL;
	if(DistanceType==0)
	{
		return DistanceStandard(IV1,IV2);
	}
	else if(DistanceType==1||DistanceType==2){
		return DistanceAdjusted(IV1,IV2,DistanceParameter);
	}
	else
	{
		cout << "!!ERROR!! !!!Distance!!!"
		     << " DistanceType defined by the user isn't recognized by the code"<<endl;

		exit(1);
	}	
DBGL;
}


double RelativDistance(IsotopicVector IV1, IsotopicVector IV2 )
{
DBGL;
	double d2 = 0;
	
	IsotopicVector IVtmp = IV1 + IV2;
	map<ZAI ,double> IVtmpIsotopicQuantity = IVtmp.GetIsotopicQuantity();
	
	double Z1total = 0;
	double Z2total = 0;
	map<ZAI ,double >::iterator it;
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		Z1total += IV1.GetZAIIsotopicQuantity( (*it).first );
		Z2total += IV2.GetZAIIsotopicQuantity( (*it).first );
	}	
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		double Z1 = IV1.GetZAIIsotopicQuantity( (*it).first );
		double Z2 = IV2.GetZAIIsotopicQuantity( (*it).first );
		d2 += pow( (Z1/Z1total - Z2/Z2total) , 2 );
	}
	
DBGL;
	return sqrt(d2);
}


IsotopicVector operator+(IsotopicVector const& IVa, IsotopicVector const& IVb)
{
DBGL;
	IsotopicVector IVtmp;
	IVtmp = IVa;
	return IVtmp += IVb;
DBGL;
}

//________________________________________________________________________
IsotopicVector operator-(IsotopicVector const& IVa, IsotopicVector const& IVb)
{
DBGL;
	IsotopicVector IVtmp;
	IVtmp = IVa;
	return IVtmp -= IVb;
DBGL;
}


//________________________________________________________________________
IsotopicVector operator*(ZAI const& zai, double F)
{
DBGL;
	IsotopicVector IVtmp;
		
	IVtmp.Add( zai, F);
	return IVtmp;
DBGL;
}


//________________________________________________________________________
IsotopicVector operator/(ZAI const& zai, double F)
{
	IsotopicVector IVtmp;
		
	IVtmp.Add( zai, 1./F);
		
	return IVtmp;
}


//________________________________________________________________________
IsotopicVector operator*(double F, IsotopicVector const& IVA) {return IVA*F;}

//________________________________________________________________________
IsotopicVector operator*(IsotopicVector const& IVA, double F)
{
DBGL;
	IsotopicVector IV = IVA;
	IV.Multiply(F);
	return IV;
DBGL;
}

IsotopicVector operator*(double F, ZAI const& zai) {return zai*F;}

//________________________________________________________________________
IsotopicVector operator/(IsotopicVector const& IVA, double F)
{
DBGL;
	IsotopicVector IV = IVA;
	IV.Multiply(1./F);
	return IV;
DBGL;
}


//____________________________InClass Operator____________________________

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator+=(const IsotopicVector& IVa)
{
DBGL;
	Add(IVa);
	return *this;
DBGL;
}

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator-=(const IsotopicVector& IVa)
{
DBGL;
	Remove(IVa);
	return *this;
DBGL;
}

bool IsotopicVector::operator<(const IsotopicVector& isotopicvector) const
{
DBGL;
	if( Norme(*this) != Norme(isotopicvector) )
		return Norme(*this) < Norme(isotopicvector);
	else if( (*this).GetIsotopicQuantity().size() != isotopicvector.GetIsotopicQuantity().size() )
		return (*this).GetIsotopicQuantity().size() < isotopicvector.GetIsotopicQuantity().size();
	else
	{
		map<ZAI ,double>::iterator it;
		map<ZAI ,double>::iterator it2 = isotopicvector.GetIsotopicQuantity().begin();
		map<ZAI ,double> IsotopicQuantity = (*this).GetIsotopicQuantity();
		for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++ )
		{
			if( (*it).first != (*it2).first )
				return (*it).first < (*it2).first;
			else it2++;
		}
		return false;
	}
DBGL;
}


//________________________________________________________________________
//________________________Constructor & Destructor________________________
//________________________________________________________________________
IsotopicVector::IsotopicVector()
{
DBGL;
DBGL;
}


//________________________________________________________________________
IsotopicVector::~IsotopicVector()
{
DBGL;
}



//________________________________________________________________________
//_____________________________General Method_____________________________
//________________________________________________________________________
void IsotopicVector::Clear()
{
DBGL;
	fIsotopicQuantityNeeded.clear();
	fIsotopicQuantity.clear();
DBGL;
}
//________________________________________________________________________
void IsotopicVector::ClearNeed()
{
DBGL;
	fIsotopicQuantityNeeded.clear();
DBGL;
}

//________________________________________________________________________
void IsotopicVector::Multiply(double factor)
{
DBGL;
	map<ZAI ,double >::iterator it;
	for( it = fIsotopicQuantity.begin(); it != fIsotopicQuantity.end(); it++)
		(*it).second = (*it).second * factor;
DBGL;
	
}

//________________________________________________________________________
void IsotopicVector::Add(const ZAI& zai, double quantity)
{
DBGL;
	if( ceil(quantity*1e25) - quantity*1e25 >  quantity*1e25 - floor(quantity*1e25) )
		quantity = floor(quantity*1e25)*1/1e25;
	else	quantity = ceil(quantity*1e25)*1/1e25;
	

	if(quantity > 0)
	{
		pair<map<ZAI, double>::iterator, bool> IResult;
		IResult = fIsotopicQuantity.insert( pair<ZAI ,double>(zai, quantity));
		if(IResult.second == false)
			IResult.first->second += quantity;
	}

DBGL;
}
//________________________________________________________________________

void IsotopicVector::Add(const IsotopicVector& isotopicvector)
{
DBGL;
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Add( (*it).first, (*it).second);
DBGL;

}
//________________________________________________________________________

void IsotopicVector::Add(const map<ZAI ,double>& quantity)
{
DBGL;
	map<ZAI ,double> isotopicquantity = quantity;
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Add( (*it).first, (*it).second);
DBGL;

}


//________________________________________________________________________
void IsotopicVector::Remove(const ZAI& zai, double quantity)
{
DBGL;

	map<ZAI ,double>::iterator it;
	it = fIsotopicQuantity.find(zai);

	if( ceil(quantity*1e25) - quantity*1e25 >  quantity*1e25 - floor(quantity*1e25) )
		quantity = floor(quantity*1e25)*1/1e25;
	else	quantity = ceil(quantity*1e25)*1/1e25;
	
	if(quantity > 0)
	{	
		if ( it != fIsotopicQuantity.end() ) 
		{
			if (it->second > quantity) 
				it->second = it->second - quantity;
			else 
			{
				Need(zai, quantity - it->second );
				it->second = 0;
			}
		}
		else
		{
			Need(zai, quantity);
		}
	}
DBGL;

}

//________________________________________________________________________
void IsotopicVector::Remove(const IsotopicVector& isotopicvector)
{
DBGL;
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Remove( (*it).first, (*it).second);
DBGL;
}

//________________________________________________________________________
void IsotopicVector::Need(const ZAI& zai, double quantity)
{
DBGL;
	if( ceil(quantity*1e25) - quantity*1e25 >  quantity*1e25 - floor(quantity*1e25) )
		quantity = floor(quantity*1e25)*1/1e25;
	else	quantity = ceil(quantity*1e25)*1/1e25;
	

	pair<map<ZAI, double>::iterator, bool> IResult;
	if(quantity > 0)
	{
		IResult = fIsotopicQuantityNeeded.insert( pair<ZAI ,double>(zai, quantity));
		if(IResult.second == false)
			IResult.first->second += quantity;
	}
DBGL;
	

}

//________________________________________________________________________
void IsotopicVector::Need(const IsotopicVector& isotopicvector)
{
DBGL;
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Need( (*it).first, (*it).second);
DBGL;
}


//________________________________________________________________________
double	IsotopicVector::GetZAIIsotopicQuantity(const ZAI& zai) const
{
DBGL;
	map<ZAI ,double> IsotopicQuantity = fIsotopicQuantity;
	
	map<ZAI ,double>::iterator it;
	it = IsotopicQuantity.find(zai);
	DBGL;
	
	if ( it != IsotopicQuantity.end() ) 
	{
		return it->second;
	}
	else
	{
		return 0;
	}
}

//________________________________________________________________________
double	IsotopicVector::GetZAIIsotopicQuantity(const int z, const int a, const int i) const
{
DBGL;
	ZAI zai(z, a, i);
	return GetZAIIsotopicQuantity(zai);
}

IsotopicVector	IsotopicVector::GetSpeciesComposition(int z) const
{
DBGL;
	IsotopicVector IV;
	map<ZAI ,double > IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)

		if( (*it).first.Z() == z  )  
			IV += (*it).first * (*it).second;

	return IV;
DBGL;
}

IsotopicVector	IsotopicVector::GetActinidesComposition() const
{
	DBGL;
	IsotopicVector IV;
	for (int i = 0; i <12; i++)
		IV += GetSpeciesComposition(89+i);
	return IV;
	DBGL;
}

vector<int> IsotopicVector::GetChemicalSpecies() const
{
DBGL;
	vector<int> ChemicalSpecies;
	
	map<ZAI ,double > IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
		if( (int)ChemicalSpecies.size() ==0 || (*it).first.Z() != ChemicalSpecies.back() )  
			ChemicalSpecies.push_back((*it).first.Z());

DBGL;
	return ChemicalSpecies;
}


//________________________________________________________________________
void IsotopicVector::Write(string filename, cSecond time) const
{
	ofstream IVfile(filename.c_str(), ios_base::app);		// Open the File
	if(!IVfile)
	cout << "!!Warning!! !!!IsotopicVector!!! \n Can't open \"" << filename << "\"\n" << endl;
	
	if(time != -1)
		IVfile << "Time "<< time/365.25/3600./24. << endl;

	map<ZAI ,double> IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for(it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
	{
		IVfile << (*it).first.Z() << " ";
		IVfile << (*it).first.A() << " ";
		IVfile << (*it).first.I() << " ";
		IVfile << (*it).second << " " << endl;
	}
	IVfile << endl;
}
//________________________________________________________________________
void IsotopicVector::Print(string option) const 
{
DBGL;
	cout << "**************************" << endl;
	cout << "*Isotopic Vector Property*" << endl;
	cout << "**************************" << endl << endl;

	bool QuantityPrint = false;
	bool DBPrint = false;

	QuantityPrint = true;
	
	if(QuantityPrint)
	{
		cout << "*Isotopic Vector Quantity*" << endl;
		map<ZAI ,double> IsotopicQuantity = GetIsotopicQuantity();
		map<ZAI ,double >::iterator it;
		for(it = IsotopicQuantity.begin();it != IsotopicQuantity.end(); it++)
		{
			cout << (*it).first.Z() << " ";
			cout << (*it).first.A() << " ";
			cout << (*it).first.I() << " ";
			cout << ": " << (*it).second;
			cout << endl;
		}
		cout << endl;
		cout << "*Isotopic Vector Quantity Needed*" << endl;
		map<ZAI ,double> IsotopicQuantityNeeded = GetIsotopicQuantityNeeded();
		for(it = IsotopicQuantityNeeded.begin(); it != IsotopicQuantityNeeded.end(); it++)
		{
			cout << (*it).first.Z() << " ";
			cout << (*it).first.A() << " ";
			cout << (*it).first.I() << " ";
			cout << ": " << (*it).second;
			cout << endl;
		}
		cout << endl;
	}
	if(DBPrint)
	{
		cout << "****Isotopic Vector DB****" << endl;
	}
DBGL;
}






