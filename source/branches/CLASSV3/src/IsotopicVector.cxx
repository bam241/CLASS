#include "IsotopicVector.hxx"

#include "LogFile.hxx"


#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "CLASSHeaders.hxx"
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

	IsotopicVector IV;

	return Distance(IV1, IV, DistanceType,DistanceParameter);

}
double DistanceStandard(IsotopicVector IV1, IsotopicVector IV2)
{

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

}
double DistanceAdjusted(IsotopicVector IV1, IsotopicVector IV2, IsotopicVector DistanceParameter)
{

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

}



double Distance(IsotopicVector IV1, IsotopicVector IV2 ,int DistanceType, IsotopicVector DistanceParameter)
{

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

}


double RelativDistance(IsotopicVector IV1, IsotopicVector IV2 )
{

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


	return sqrt(d2);
}


IsotopicVector operator+(IsotopicVector const& IVa, IsotopicVector const& IVb)
{

	IsotopicVector IVtmp;
	IVtmp = IVa;
	return IVtmp += IVb;
}

//________________________________________________________________________
IsotopicVector operator-(IsotopicVector const& IVa, IsotopicVector const& IVb)
{

	IsotopicVector IVtmp;
	IVtmp = IVa;
	return IVtmp -= IVb;
}


//________________________________________________________________________
IsotopicVector operator*(ZAI const& zai, double F)
{

	IsotopicVector IVtmp;

	IVtmp.Add( zai, F);
	return IVtmp;
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

	IsotopicVector IV = IVA;
	IV.Multiply(F);
	return IV;
}

//________________________________________________________________________
IsotopicVector operator*(double F, ZAI const& zai)
{
	return zai*F;
}

//________________________________________________________________________
IsotopicVector operator*(IsotopicVector const& IVa, IsotopicVector const& IVb)
{

	IsotopicVector IVtmp;
	IVtmp = IVa;
	return IVtmp *= IVb;
}

//________________________________________________________________________
IsotopicVector operator/(IsotopicVector const& IVA, double F)
{

	IsotopicVector IV = IVA;
	IV.Multiply(1./F);
	return IV;
}


//____________________________InClass Operator____________________________

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator+=(const IsotopicVector& IVa)
{

	Add(IVa);
	return *this;

}

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator-=(const IsotopicVector& IVa)
{

	Remove(IVa);
	return *this;

}
//________________________________________________________________________
IsotopicVector& IsotopicVector::operator*=(const IsotopicVector& IVa)
{
	map<ZAI, double> IVA_isotopicquantity = IVa.GetIsotopicQuantity();

	map<ZAI, double> isotopicquantity = (*this).GetIsotopicQuantity(); // get the isotopic quantity to loop on it
	map<ZAI, double>::iterator isotopicIT;				  // iterator on a isotopic quantity map

	for(isotopicIT = isotopicquantity.begin(); isotopicIT != isotopicquantity.end(); isotopicIT++) // loop on the isotopicquantity...
	{
		map<ZAI, double>::iterator IVa_isotopicIT = IVA_isotopicquantity.find( (*isotopicIT).first );

		if(IVa_isotopicIT != IVA_isotopicquantity.end())
		(*isotopicIT).second *= (*IVa_isotopicIT).second;
		else
		(*isotopicIT).second *= 0;

	}

	return *this;

}

//________________________________________________________________________
bool IsotopicVector::operator<(const IsotopicVector& isotopicvector) const
{

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

}


//________________________________________________________________________
//________________________Constructor & Destructor________________________
//________________________________________________________________________
IsotopicVector::IsotopicVector()
{


}


//_____________________________________________________GetSpeciesComposition___________________
IsotopicVector::~IsotopicVector()
{
	fIsotopicQuantity.clear();
	fIsotopicQuantityNeeded.clear();
}



//________________________________________________________________________
//_____________________________General Method_____________________________
//________________________________________________________________________
void IsotopicVector::Clear()
{

	fIsotopicQuantityNeeded.clear();
	fIsotopicQuantity.clear();

}
//________________________________________________________________________
void IsotopicVector::ClearNeed()
{

	fIsotopicQuantityNeeded.clear();

}

//________________________________________________________________________
void IsotopicVector::Multiply(double factor)
{

	map<ZAI ,double >::iterator it;
	for( it = fIsotopicQuantity.begin(); it != fIsotopicQuantity.end(); it++)
	(*it).second = (*it).second * factor;
	for( it = fIsotopicQuantityNeeded.begin(); it != fIsotopicQuantityNeeded.end(); it++)
	(*it).second = (*it).second * factor;


}

//________________________________________________________________________

double IsotopicVector::GetSumOfAll() const
{
	double Sum = 0;
	map<ZAI ,double >::iterator it;
	map<ZAI ,double > isotopicquantity = GetIsotopicQuantity();
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	Sum += (*it).second;

	return Sum;

}
//________________________________________________________________________
void IsotopicVector::Add(const ZAI& zai, double quantity)
{

	if( ceil(quantity*1e25) - quantity*1e25 >  quantity*1e25 - floor(quantity*1e25) )
	quantity = floor(quantity*1e25)*1/1e25;
	else	quantity = ceil(quantity*1e25)*1/1e25;


	if(quantity > 0)
	{
		pair<map<ZAI, double>::iterator, bool> IResult;
		IResult = fIsotopicQuantity.insert( pair<ZAI ,double>(zai, quantity));
		if(!IResult.second)
		IResult.first->second += quantity;
	}


}
//________________________________________________________________________

void IsotopicVector::Add(const IsotopicVector& isotopicvector)
{

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	Add( (*it).first, (*it).second);


}
//________________________________________________________________________

void IsotopicVector::Add(const map<ZAI ,double>& quantity)
{

	map<ZAI ,double> isotopicquantity = quantity;
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	Add( (*it).first, (*it).second);


}


//________________________________________________________________________
void IsotopicVector::Remove(const ZAI& zai, double quantity)
{


	map<ZAI ,double>::iterator it;
	it = fIsotopicQuantity.find(zai);

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


}

//________________________________________________________________________
void IsotopicVector::Remove(const IsotopicVector& isotopicvector)
{

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	Remove( (*it).first, (*it).second);

}

//________________________________________________________________________
void IsotopicVector::Need(const ZAI& zai, double quantity)
{
	pair<map<ZAI, double>::iterator, bool> IResult;
	if(quantity > 0)
	{
		IResult = fIsotopicQuantityNeeded.insert( pair<ZAI ,double>(zai, quantity));
		if(!IResult.second)
		IResult.first->second += quantity;
	}


}

//________________________________________________________________________
void IsotopicVector::Need(const IsotopicVector& isotopicvector)
{

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	Need( (*it).first, (*it).second);

}


//________________________________________________________________________
double	IsotopicVector::GetZAIIsotopicQuantity(const ZAI& zai) const
{

	map<ZAI ,double> IsotopicQuantity = fIsotopicQuantity;

	map<ZAI ,double>::iterator it;
	it = IsotopicQuantity.find(zai);


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

	ZAI zai(z, a, i);
	return GetZAIIsotopicQuantity(zai);
}

IsotopicVector	IsotopicVector::GetSpeciesComposition(int z) const
{

	IsotopicVector IV;
	map<ZAI ,double > IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
		if( (*it).first.Z() == z  )
			IV += (*it).first * (*it).second;

	return IV;

}


IsotopicVector	IsotopicVector::GetThisComposition(IsotopicVector IV) const
{
	IsotopicVector IVtmp;
	map<ZAI ,double > IsotopicQuantity = IV.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
	{
		map<ZAI ,double >::iterator it2 = GetIsotopicQuantity().find((*it).first);
		if(it2 != GetIsotopicQuantity().end())
			IVtmp += (*it2).first * (*it2).second ;

	}

	return IVtmp;
	
}
//________________________________________________________________________
double IsotopicVector::GetTotalMass() const
{
	double TotalMass = 0;
	double AVOGADRO=6.02214129e23;
	map<ZAI ,double >::iterator it;
	map<ZAI ,double > isotopicquantity = GetIsotopicQuantity();
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		ZAI zai = ((*it).first);
		double MolarMass = zai.GetMass();
		TotalMass += (*it).second/AVOGADRO * MolarMass ;
	}



	return TotalMass*1e-6;//in tons

}
//________________________________________________________________________

double IsotopicVector::MeanMolar() const
{
	double MeanMolar = 0;
	map<ZAI ,double >::iterator it;
	map<ZAI ,double > isotopicquantity = GetIsotopicQuantity();

	double NTot=0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	NTot+= (*it).second ;

	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		ZAI zai = ((*it).first);
		double MolarMass = zai.GetMass();
		MeanMolar += (*it).second/NTot * MolarMass ;
	}


	return MeanMolar;

}
//________________________________________________________________________

vector<ZAI> IsotopicVector::GetZAIList() const
{

	map<ZAI ,double > IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	vector<ZAI> zailist;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
	zailist.push_back( (*it).first );

	return zailist;

}

IsotopicVector	IsotopicVector::GetActinidesComposition() const
{

	IsotopicVector IV;
	for (int i = 0; i <12; i++)
	IV += GetSpeciesComposition(89+i);
	return IV;

}

vector<int> IsotopicVector::GetChemicalSpecies() const
{

	vector<int> ChemicalSpecies;

	map<ZAI ,double > IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
	if( (int)ChemicalSpecies.size() ==0 || (*it).first.Z() != ChemicalSpecies.back() )
	ChemicalSpecies.push_back((*it).first.Z());


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

}


void IsotopicVector::PrintList(string option) const
{
	bool QuantityPrint = false;
	bool DBPrint = false;

	QuantityPrint = true;

	if(QuantityPrint)
	{
		map<ZAI ,double> IsotopicQuantity = GetIsotopicQuantity();
		map<ZAI ,double >::iterator it;
		for(it = IsotopicQuantity.begin();it != IsotopicQuantity.end(); it++)
		{
			cout << (*it).first.Z() << " ";
			cout << (*it).first.A() << " ";
			cout << (*it).first.I() << " ";
			cout << endl;
		}
		cout << endl;

	}
	if(DBPrint)
	{
		cout << "****Isotopic Vector DB****" << endl;
	}

}




