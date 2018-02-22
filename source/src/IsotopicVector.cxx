#include "IsotopicVector.hxx"

#include "CLASSLogger.hxx"
#include "CLASSConstante.hxx"


#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
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

ClassImp(IsotopicVector)

//________________________________________________________________________
//________________________Constructor & Destructor________________________
//________________________________________________________________________
IsotopicVector::IsotopicVector ()
{ ; }

IsotopicVector::IsotopicVector ( IsotopicVector const& IVa ) :
	fIsotopicQuantity(IVa.fIsotopicQuantity) , fIsotopicQuantityNeeded(IVa.fIsotopicQuantityNeeded)
{ ; }

//_____________________________________________________GetSpeciesComposition___________________
IsotopicVector::~IsotopicVector ()
{
	//fIsotopicQuantity.clear();
	//fIsotopicQuantityNeeded.clear();
}


//____________________________InClass Operator____________________________

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator+= (const IsotopicVector& IVa)
{
	const_iterator end = IVa.end();
	for ( const_iterator it=IVa.begin() ; it!=end ; ++it )
		{ Add( *it ); }

	return *this;
}

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator-= (const IsotopicVector& IVa)
{
	const_iterator end = IVa.end();
	for ( const_iterator it=IVa.begin() ; it!=end ; ++it )
		{ Remove( *it ); }

	return *this;
}
//________________________________________________________________________
IsotopicVector& IsotopicVector::operator*= (const IsotopicVector& IVa)
{
	iterator f_end = fIsotopicQuantity.end();

	const_iterator IVa_end = IVa.end();
	const_iterator jt;

	for( iterator it = fIsotopicQuantity.begin() ; it != f_end ; ++it )
	{
		jt = IVa.find( it->first );

		if ( jt != IVa_end )
			{ it->second *=  jt->second; }
		else
			{ it->second = 0; }
	}

	return *this;
}

//________________________________________________________________________
IsotopicVector& IsotopicVector::operator*= (const double& factor)
{
	Multiply(factor);
	
	return *this;
}

//________________________________________________________________________
bool IsotopicVector::operator < ( const IsotopicVector& isotopicvector ) const
{
	if ( Norme(*this) != Norme(isotopicvector) )
		{ return Norme(*this) < Norme(isotopicvector); }
	else if ( (*this).GetIsotopicQuantity().size() != isotopicvector.GetIsotopicQuantity().size() )
		{ return (*this).GetIsotopicQuantity().size() < isotopicvector.GetIsotopicQuantity().size(); }
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
	iterator it;

	for( it = fIsotopicQuantity.begin(); it != fIsotopicQuantity.end(); it++)
		{ (*it).second = (*it).second * factor; }

	for( it = fIsotopicQuantityNeeded.begin(); it != fIsotopicQuantityNeeded.end(); it++)
		{ (*it).second = (*it).second * factor; }
}

//________________________________________________________________________
double IsotopicVector::GetSumOfAll () const
{
	return std::accumulate(
			fIsotopicQuantity.begin() , fIsotopicQuantity.end() , // input iterator
			0.0 ,                         // first value for the sum
			[]( const double sum , const std::pair<ZAI,double> & zaiQ ) { return zaiQ.second + sum; }
		);
}

//________________________________________________________________________
void IsotopicVector::Add(const ZAI& zai, double q)
{
	if( ceil(q*1e50) - q*1e50 >  q*1e50 - floor(q*1e50) )
		{ q = floor(q*1e50)*1/1e50; }
	else
		{ q = ceil(q*1e50)*1/1e50; }

	if ( q > 0 )
		{ fIsotopicQuantity[zai] += q; }
}

//________________________________________________________________________
void IsotopicVector::Add(const IsotopicVector& isotopicvector)
{
	Add( isotopicvector.fIsotopicQuantity );
}
//________________________________________________________________________
void IsotopicVector::Add(const map<ZAI ,double>& quantity)
{
	const_iterator end = quantity.end();
	for ( const_iterator it = quantity.begin() ; it!=end ; ++it )
	{
		Add( it->first , it->second );
	}
}


//________________________________________________________________________
void IsotopicVector::Remove(const ZAI& zai, double quantity)
{
	iterator it = fIsotopicQuantity.find(zai);

	if ( quantity > 0 )
	{
		if ( it != fIsotopicQuantity.end() )
		{
			if ( it->second > quantity )
				{ it->second = it->second - quantity; }
			else
			{
				if( (it->second - quantity)/it->second > 1e-16 ) // to fit with double precision : 16 digits
					{ Need( zai , quantity - it->second ); }
				it->second = 0;
			}
		}
		else
		{
			Need(zai, quantity);
		}
	}

	if ( it->second ==  0 )
		{ fIsotopicQuantity.erase(it); }
}

//________________________________________________________________________
void IsotopicVector::Remove(const IsotopicVector& IVa)
{
	const_iterator IVa_end = IVa.end();
	
	iterator f_end = fIsotopicQuantity.end();
	iterator jt;
	for ( const_iterator it=IVa.begin() ; it!=IVa_end ; ++it )
	{
		jt = fIsotopicQuantity.find( it->first );
		if ( jt != f_end )
		{
			jt->second -= it->second;
			if ( jt->second <= 0 )
			{
				fIsotopicQuantity.erase(jt);
				f_end = fIsotopicQuantity.end();
			}
		}
	}
}

//________________________________________________________________________

void IsotopicVector::ApplyZAIThreshold(int z)
{
	double quantity = 0.0;
	iterator it = fIsotopicQuantity.begin();

	while ( it != fIsotopicQuantity.end() )
	{
		if ( it->first.Z() < z )
		{
			quantity += it->second;
			it = fIsotopicQuantity.erase(it++); //becarefull of postcrement incrementation
		}
		else
		{
			++it;
		}
	}
	fIsotopicQuantity.insert( std::pair<ZAI,double>( ZAI(-2,-2,-2) , quantity ) );
}




//________________________________________________________________________
void IsotopicVector::Need(const ZAI& zai, double quantity)
{
	if(quantity < 0.5) quantity = 0;

	pair<map<ZAI, double>::iterator, bool> IResult;
	if(quantity > 0)
	{	cout << "Negative quantity : " << quantity  << " for ZAI " << zai.Z() << " " << zai.A() << " " << zai.I() << " in this IsotopicVector" << endl;
		exit(0);
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
	const_iterator it = fIsotopicQuantity.find(zai);
	if ( it != fIsotopicQuantity.end() )
		{ return it->second; }

	return 0;
}

//________________________________________________________________________
double	IsotopicVector::GetZAIIsotopicQuantity(const int z, const int a, const int i) const
	{ return GetZAIIsotopicQuantity(ZAI(z,a,i)); }

//________________________________________________________________________
IsotopicVector	IsotopicVector::GetSpeciesComposition ( const int z ) const
{
	IsotopicVector tmp;
	const_iterator end = fIsotopicQuantity.end();

	for ( const_iterator it=fIsotopicQuantity.begin() ; it!=end ; ++it )
	{
		if ( it->first.Z() == z )
		{
			tmp.fIsotopicQuantity.insert( *it );
		}
	}

	return tmp;
}
//________________________________________________________________________
IsotopicVector	IsotopicVector::GetThisChemicalComposition(IsotopicVector IVa) const
{
	IsotopicVector tmp;
	vector<int> ZList;
	
	ZList = IVa.GetChemicalSpecies();
	
	for ( int i = 0 ; i<ZList.size() ; i++)
	{		
		for( const_iterator jt = fIsotopicQuantity.begin() ; jt != fIsotopicQuantity.end() ; jt ++)
		if (jt->first.Z() == ZList[i])
		{
			tmp.fIsotopicQuantity.insert(*jt);
		}
		
	}

	return tmp;
}

//________________________________________________________________________
IsotopicVector	IsotopicVector::GetThisComposition(IsotopicVector IVa) const
{
	IsotopicVector tmp;

	const_iterator IVa_end = IVa.end();
	const_iterator f_end = fIsotopicQuantity.end();

	const_iterator jt;
	for ( const_iterator it = IVa.begin() ; it != IVa_end ; ++it )
	{
		jt = fIsotopicQuantity.find( it->first );
		if ( jt != f_end )
			{ tmp.fIsotopicQuantity.insert( *jt ); }
	}

	return tmp;
}
//________________________________________________________________________
double IsotopicVector::GetTotalMass() const
{
	return cZAIMass.GetMass(*this);//in tons
}
//________________________________________________________________________

double IsotopicVector::GetMeanMolarMass() const
{
	return GetTotalMass() * 1e6 * AVOGADRO / GetActinidesComposition().GetSumOfAll();
}

//________________________________________________________________________
vector<ZAI> IsotopicVector::GetZAIList() const
{
	std::vector< ZAI > tmp( fIsotopicQuantity.size() );

	std::transform (
			fIsotopicQuantity.begin() , fIsotopicQuantity.end() , // input iterator
			tmp.begin()                 , // output iterator
			[]( const std::pair<ZAI,double> & zaiQ ) { return zaiQ.first; }
		);

	return tmp;
}
//________________________________________________________________________
void IsotopicVector::Initiatlize(double val)
{
	if ( val == 0 ) { fIsotopicQuantity.clear(); }
	else
	{	
		iterator end = fIsotopicQuantity.end();
		for ( iterator it=fIsotopicQuantity.begin() ; it!=end ; ++it )
			{ it->second = val; }
	}
}

//________________________________________________________________________
IsotopicVector	IsotopicVector::GetActinidesComposition() const
{
	IsotopicVector tmp;
	const_iterator end = fIsotopicQuantity.end();

	for ( const_iterator it=fIsotopicQuantity.begin() ; it!=end ; ++it )
	{
		if ( it->first.Z() >= 89 && it->first.Z() <= 103 )
		{
			tmp.fIsotopicQuantity.insert( *it );
		}
	}

	return tmp;
}

vector<int> IsotopicVector::GetChemicalSpecies() const
{
	std::vector< int > tmp( fIsotopicQuantity.size() );

	std::transform (
			fIsotopicQuantity.begin() , fIsotopicQuantity.end() , // input iterator
			tmp.begin()                 , // output iterator
			[]( const std::pair<ZAI,double> & zaiQ ) { return zaiQ.first.Z(); }
		);

	std::vector<int>::iterator end = std::unique( tmp.begin() , tmp.end() );
	tmp.resize( end - tmp.begin() );

	return tmp;
}


//________________________________________________________________________
void IsotopicVector::Write(string filename, cSecond time) const
{
	ofstream IVfile(filename.c_str(), ios_base::app);		// Open the File
	if(!IVfile)
	cout << "!!Warning!! !!!IsotopicVector!!! \n Can't open \"" << filename << "\"\n" << endl;

	if(time != -1)
	IVfile << "Time " <<  time/cYear << endl;

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

	cout << sPrint();

}
//________________________________________________________________________
string IsotopicVector::sPrint() const
{
	stringstream ss;
	ss << "**************************" << endl;
	ss << "*Isotopic Vector Property*" << endl;
	ss << "**************************" << endl << endl;

	bool QuantityPrint = false;
	bool DBPrint = false;

	QuantityPrint = true;

	if(QuantityPrint)
	{
		ss << "*Isotopic Vector Quantity*" << endl;
		map<ZAI ,double> IsotopicQuantity = GetIsotopicQuantity();
		map<ZAI ,double >::iterator it;
		for(it = IsotopicQuantity.begin();it != IsotopicQuantity.end(); it++)
		{
			ss << (*it).first.Z() << " ";
			ss << (*it).first.A() << " ";
			ss << (*it).first.I() << " ";
			ss << ": " << (*it).second;
			ss << endl;
		}
		ss << endl;
		ss << "*Isotopic Vector Quantity Needed*" << endl;
		map<ZAI ,double> IsotopicQuantityNeeded = GetIsotopicQuantityNeeded();
		for(it = IsotopicQuantityNeeded.begin(); it != IsotopicQuantityNeeded.end(); it++)
		{
			ss << (*it).first.Z() << " ";
			ss << (*it).first.A() << " ";
			ss << (*it).first.I() << " ";
			ss << ": " << (*it).second;
			ss << endl;
		}
		ss << endl;
	}
	if(DBPrint)
	{
		ss << "****Isotopic Vector DB****" << endl;
	}
return ss.str();
}
//________________________________________________________________________
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

//________________________________________________________________________
//__________________________Operator Overlaoding__________________________
//________________________________________________________________________


IsotopicVector operator+(IsotopicVector const& IVa, IsotopicVector const& IVb)
{
	IsotopicVector IVtmp = IVa;
	return IVtmp += IVb;
}

//________________________________________________________________________
IsotopicVector operator-(IsotopicVector const& IVa, IsotopicVector const& IVb)
{

	IsotopicVector IVtmp = IVa;
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
IsotopicVector operator/(IsotopicVector const& IVA, IsotopicVector const& IVB)
{
	IsotopicVector IVAFromB  = IVA.GetThisComposition(IVB);
	unsigned int IVAFromB_Size = IVAFromB.GetZAIList().size();
	unsigned int IVA_Size      = IVA.GetZAIList().size();
	unsigned int IVB_Size      = IVB.GetZAIList().size();

	if( IVAFromB_Size < IVA_Size)
	{
		cout << "Something try to divide by zero" << endl;
		cout << "IVA / IVB All ZAI in IVA have to be present in IVB  " << endl;
		exit(0);
	}
	else if ( IVB_Size > IVAFromB_Size )
		cout << " IVA / IVB : Size of IVB is bigger than IVA. Non commun nuclei will not be present in IVresult" << endl;
	
	IsotopicVector IVresult;
	map<ZAI,double>::const_iterator end   = IVAFromB.end();
	map<ZAI,double>::const_iterator begin = IVAFromB.begin();
	map<ZAI,double>::const_iterator it;

	for ( it=begin ; it!=end ; ++it )
		IVresult.Add(it->first,it->second/IVB.GetIsotopicQuantity()[it->first]);

	return IVresult;

}

//________________________________________________________________________
IsotopicVector operator*(double F, IsotopicVector const& IVa)
	{return IVa*F;}

//________________________________________________________________________
IsotopicVector operator*(IsotopicVector const& IVa, double F)
{
	IsotopicVector IV = IVa;
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
	IVtmp *=  IVb;
	return IVtmp;
}

//________________________________________________________________________
IsotopicVector operator/(IsotopicVector const& IVa, double F)
{

	IsotopicVector IV = IVa;
	IV.Multiply(1./F);
	return IV;
}


//____________________________General Operator____________________________
//________________________________________________________________________
double 	RelativDistance ( const IsotopicVector & a, const IsotopicVector & b )
{
	double d2 = 0;

	IsotopicVector tmp = a + b;

	double Z1total = 0;
	double Z2total = 0;
	for( IsotopicVector::iterator it = tmp.begin(); it != tmp.end(); ++it )
	{
		Z1total += a.GetZAIIsotopicQuantity( it->first );
		Z2total += b.GetZAIIsotopicQuantity( it->first );
	}
	double Z1, Z2;
	for( IsotopicVector::iterator it = tmp.begin(); it != tmp.end(); ++it )
	{
		Z1 = a.GetZAIIsotopicQuantity( it->first );
		Z2 = b.GetZAIIsotopicQuantity( it->first );
		d2 += (Z1/Z1total - Z2/Z2total)*(Z1/Z1total - Z2/Z2total);
	}

	return std::sqrt(d2);
}
//____________________________________________________________________________
double 	Distance ( const IsotopicVector & a , const IsotopicVector & b , int DistanceType , const IsotopicVector & DistanceParameter )
{
	if(DistanceType == 0)
	{
		return DistanceStandard(a,b);
	}
	else if(DistanceType == 1||DistanceType == 2){
		return DistanceAdjusted(a,b,DistanceParameter);
	}
	else
	{
		std::cout << " DistanceType defined by the user isn't recognized by the code" << std::endl;

		exit(1);
	}
}
//____________________________________________________________________________
double	DistanceStandard( const IsotopicVector & a , const IsotopicVector & b )
{
	double d2 = 0.0;
	IsotopicVector tmp = a + b;
	double Z1,Z2;
	for( IsotopicVector::iterator it = tmp.begin(); it != tmp.end(); ++it )
	{
		Z1 = a.GetZAIIsotopicQuantity( it->first );
		Z2 = b.GetZAIIsotopicQuantity( it->first );
		d2 += (Z1-Z2)*(Z1-Z2);
	}
	return std::sqrt(d2);
}
//____________________________________________________________________________
double	DistanceAdjusted( const IsotopicVector & a , const IsotopicVector & b , const IsotopicVector & DistanceParameter )
{
	double d2 = 0;
	IsotopicVector tmp = a + b;
	double Z1,Z2;
	double lambda;
	for( IsotopicVector::iterator it = tmp.begin(); it != tmp.end(); ++it )
	{
		Z1 = a.GetZAIIsotopicQuantity( it->first );
		Z2 = b.GetZAIIsotopicQuantity( it->first );

		lambda = DistanceParameter.GetZAIIsotopicQuantity( it->first );
		d2 += lambda*std::abs(Z1-Z2);
	}
	return d2;
}
//____________________________________________________________________________
double 	Norme( const IsotopicVector & a , int DistanceType , const IsotopicVector & DistanceParameter )
{
	IsotopicVector zero;

	return Distance(a, zero, DistanceType,DistanceParameter);
}


