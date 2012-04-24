#include "IsotopicVector.hxx"
#include "ZAI.hxx"
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

double 	Norme(IsotopicVector IV1)
{
	DBGL;
	IsotopicVector IV; 
	return Distance(IV1, IV);
	DBGL;
}

double Distance(IsotopicVector IV1, IsotopicVector IV2 )
{
	DBGL;
	double d2 = 0;
	
	IsotopicVector IVtmp = IV1 + IV2;
	map<ZAI ,double> IVtmpIsotopicQuantity = IVtmp.GetIsotopicQuantity();
	
	for(map<ZAI ,double >::iterator it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		double Z1 = IV1.GetZAIIsotopicQuantity((*it).first);
		double Z2 = IV2.GetZAIIsotopicQuantity((*it).first);
		d2 += pow(Z1-Z2 , 2 );
	}
	
	DBGL;
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
IsotopicVector operator*(IsotopicVector const& IVA, double F)
{
	DBGL;
	IsotopicVector IV = IVA;
	IV.Multiply(F);
	return IV;
}

//________________________________________________________________________
IsotopicVector operator/(IsotopicVector const& IVA, double F)
{
	DBGL;
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
	for(map<ZAI ,double >::iterator it = fIsotopicQuantity.begin(); it != fIsotopicQuantity.end(); it++)
		(*it).second = (*it).second * factor;
	DBGL;
	
}

//________________________________________________________________________
void IsotopicVector::Add(const ZAI& zai, double quantity)
{
	DBGL;
	if( ceil(quantity*1e6) - quantity*1e6 >  quantity*1e6 - floor(quantity*1e6) )
		quantity = floor(quantity*1e6)*1/1e6;
	else	quantity = ceil(quantity*1e6)*1/1e6;
	
	pair<map<ZAI, double>::iterator, bool> IResult;
	if(quantity > 0)
	{
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
	
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Add( (*it).first, (*it).second);
	DBGL;

}
//________________________________________________________________________

void IsotopicVector::Add(const map<ZAI ,double>& quantity)
{
	DBGL;
	map<ZAI ,double> isotopicquantity = quantity;
	
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Add( (*it).first, (*it).second);
	DBGL;

}


//________________________________________________________________________
void IsotopicVector::Remove(const ZAI& zai, double quantity)
{
	DBGL;

	map<ZAI ,double>::iterator it;
	it = fIsotopicQuantity.find(zai);

	if( ceil(quantity*1e6) - quantity*1e6 >  quantity*1e6 - floor(quantity*1e6) )
		quantity = floor(quantity*1e6)*1/1e6;
	else	quantity = ceil(quantity*1e6)*1/1e6;
	
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
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Remove( (*it).first, (*it).second);
	DBGL;
}

//________________________________________________________________________
void IsotopicVector::Need(const ZAI& zai, double quantity)
{
	DBGL;
	if( ceil(quantity*1e6) - quantity*1e6 >  quantity*1e6 - floor(quantity*1e6) )
		quantity = floor(quantity*1e6)*1/1e6;
	else	quantity = ceil(quantity*1e6)*1/1e6;
	

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
	
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		Need( (*it).first, (*it).second);
	DBGL;
}


//________________________________________________________________________
double	IsotopicVector::GetZAIIsotopicQuantity(const ZAI& zai) const
{
	DBGL;
	map<ZAI ,double>::iterator it;
	it = GetIsotopicQuantity().find(zai);

	DBGL;
	if ( it != GetIsotopicQuantity().end() ) 
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
	
//	map<ZAI ,double>::iterator it;
//	it = GetIsotopicQuantity().find(zai);

//	DBGL;
//	if ( it != GetIsotopicQuantity().end() ) 
//	{
//		return it->second;	
//	}	
//	else	
//	{	
//		return 0;
//	}
}

IsotopicVector	IsotopicVector::GetAtomicComposition(int z) const
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

vector<int> IsotopicVector::GetAtomicSpecies() const
{
	DBGL;
	vector<int> AtomicSpecies;
	
	map<ZAI ,double > IsotopicQuantity = GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
		if( (int)AtomicSpecies.size() ==0 || (*it).first.Z() != AtomicSpecies.back() )  
			AtomicSpecies.push_back((*it).first.Z());

	DBGL;
	return AtomicSpecies;
}


//________________________________________________________________________
void IsotopicVector::Write(string filename, long int time) const 
{
	ofstream IVfile(filename.c_str(), ios_base::app);		// Open the File
	if(!IVfile)
	cout << "!!Warning!! !!!IsotopicVector!!! \n Can't open \"" << filename << "\"\n" << endl;

	IVfile << time << " ";
	
	map<ZAI ,double> IsotopicQuantity = GetIsotopicQuantity();
	for(map<ZAI ,double >::iterator it = IsotopicQuantity.begin();
					it != IsotopicQuantity.end(); it++)
	{
		IVfile << (*it).first.Z() << " ";
		IVfile << (*it).first.A() << " ";
		IVfile << (*it).first.I() << " ";
		IVfile << (*it).second << " ";
		
		if((*it).first.Z()>70)
		{
			char Z[33];
			sprintf(Z,"%d",(*it).first.Z());
			char A[33];
			sprintf(A,"%d",(*it).first.A());
			char I[33];
			sprintf(I,"%d",(*it).first.I());
		
		
			string filenameZAI = filename+ "_DIR" +"/" + Z + "_" + A + "_"+ I;
			string cmd = " mkdir " + filename+ "_DIR" +" 2>/dev/null";
			system(cmd.c_str());
			ofstream IVfileZAI(filenameZAI.c_str(), ios_base::app);		// Open the File

			if(!IVfileZAI)
				cout << "!!Warning!! !!!IsotopicVector!!! \n Can't open \"" << filenameZAI << "\"\n" << endl;
		
			IVfileZAI << time << " " << (*it).second << endl;;
			IVfileZAI.close();		
		}
	
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
	
	for(int i = 0; i < (int)option.length(); i++)
	{
//		if(option[i] == 'd')	DBPrint = true;	
	}
	
	if(QuantityPrint)
	{
		cout << "*Isotopic Vector Quantity*" << endl;
		map<ZAI ,double> IsotopicQuantity = GetIsotopicQuantity();
		for(map<ZAI ,double >::iterator it = IsotopicQuantity.begin();
						it != IsotopicQuantity.end(); it++)
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
		for(map<ZAI ,double >::iterator it = IsotopicQuantityNeeded.begin();
		    				it != IsotopicQuantityNeeded.end(); it++)
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






