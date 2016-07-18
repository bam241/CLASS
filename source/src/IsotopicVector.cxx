#include "IsotopicVector.hxx"

#include "CLASSLogger.hxx"
#include "CLASSConstante.hxx"

ClassImp(IsotopicVector)

///// CONSTRUCTOR ////////////////////////////////////////////////////////////
IsotopicVector::IsotopicVector () 
{ ; }
//____________________________________________________________________________
IsotopicVector::IsotopicVector ( const IsotopicVector & a ) :
	fdata(a.fdata)
{ ; }

///// DESTRUCTOR /////////////////////////////////////////////////////////////
IsotopicVector::~IsotopicVector ()
{ ; }

///// OPERATOR ///////////////////////////////////////////////////////////////
IsotopicVector & IsotopicVector::operator = ( const IsotopicVector & a )
{
	if ( this != &a )
	{
		fdata = a.fdata;
		fIsotopicQuantityNeeded = a.fIsotopicQuantityNeeded;
	}
	return *this;
}

//____________________________________________________________________________
double IsotopicVector::operator [] ( const ZAI & zai ) const
{
	const_iterator it = fdata.find(zai);
	if ( it != fdata.end() )
		{ return it->second; }

	return 0;
}
//____________________________________________________________________________
double & IsotopicVector::operator [] ( const ZAI & zai )
	{ return fdata[zai]; }

//____________________________________________________________________________
IsotopicVector & IsotopicVector::operator += ( const IsotopicVector & a )
{
	const_iterator end = a.end();
	for ( const_iterator it=a.begin() ; it!=end ; ++it )
		{ Add( *it ); }

	return *this;
}
//____________________________________________________________________________
IsotopicVector & IsotopicVector::operator -= ( const IsotopicVector & a )
{
	const_iterator end = a.end();
	for ( const_iterator it=a.begin() ; it!=end ; ++it )
		{ Remove( *it ); }

	return *this;
}
//____________________________________________________________________________
IsotopicVector & IsotopicVector::operator *= ( const IsotopicVector & a )
{
	const_iterator a_end = a.end(); // just to prevent to much call of a function
	const_iterator f_end = fdata.end(); // just to prevent to much call of a function

	iterator tmp;
	for ( const_iterator it=a.begin() ; it!=a_end ; ++it )
	{
		tmp = fdata.find(it->first);
		if ( tmp != f_end )
			{ tmp->second *= it->second; }
	}

	return *this;
}
//____________________________________________________________________________
IsotopicVector & IsotopicVector::operator *= ( double f )
{
	const_iterator end = fdata.end();
	for ( iterator it=fdata.begin() ; it!=end ; ++it )
	{
		it->second *= f;
	}

	return *this;
}

///// GETTER /////////////////////////////////////////////////////////////////
std::map<ZAI,double> IsotopicVector::GetIsotopicQuantity() const
	{ return fdata; }
//____________________________________________________________________________
IsotopicVector IsotopicVector::GetSpeciesComposition ( const int a ) const
{
	IsotopicVector tmp;
	const_iterator end = fdata.end();

	for ( const_iterator it=fdata.begin() ; it!=end ; ++it )
	{
		if ( it->first.Z() == a )
		{
			tmp.fdata.insert( *it );
		}
	}

	return tmp;
}
//____________________________________________________________________________
IsotopicVector IsotopicVector::GetThisComposition ( const IsotopicVector & a ) const
{
	IsotopicVector tmp;

	const_iterator a_end = a.end();
	const_iterator f_end = fdata.end();

	const_iterator jt;
	for ( const_iterator it = a.begin() ; it != a_end ; ++it )
	{
		jt = fdata.find( it->first );
		tmp.fdata.insert( *jt );
	}

	return tmp;
}
//____________________________________________________________________________
std::vector< ZAI > IsotopicVector::GetZAIList () const
{
	std::vector< ZAI > tmp( fdata.size() );

	std::transform (
			fdata.begin() , fdata.end() , // input iterator
			tmp.begin()                 , // output iterator
			[]( const std::pair<ZAI,double> & zaiQ ) { return zaiQ.first; }
		);

	return tmp;
}
//____________________________________________________________________________
IsotopicVector IsotopicVector::GetActinidesComposition () const
{
	IsotopicVector tmp;
	const_iterator end = fdata.end();

	for ( const_iterator it=fdata.begin() ; it!=end ; ++it )
	{
		if ( it->first.Z() >= 89 && it->first.Z() <= 103 )
		{
			tmp.fdata.insert( *it );
		}
	}

	return tmp;
}

//____________________________________________________________________________
double IsotopicVector::GetZAIIsotopicQuantity ( const short int z , const short int a , const short int i ) const
{
	const_iterator it = fdata.find(ZAI(z,a,i));
	if ( it != fdata.end() )
		{ return it->second; }

	return 0;
}
//____________________________________________________________________________
double IsotopicVector::GetZAIIsotopicQuantity ( const ZAI & zai ) const
{
	const_iterator it = fdata.find(zai);
	if ( it != fdata.end() )
		{ return it->second; }

	return 0;
}
//____________________________________________________________________________
double IsotopicVector::GetQuantity ( const short int z , const short int a , const short int i ) const
{
	const_iterator it = fdata.find(ZAI(z,a,i));
	if ( it != fdata.end() )
		{ return it->second; }

	return 0;
}
//____________________________________________________________________________
double IsotopicVector::GetQuantity ( const ZAI & zai ) const
{
	const_iterator it = fdata.find(zai);
	if ( it != fdata.end() )
		{ return it->second; }

	return 0;
}

//____________________________________________________________________________
double IsotopicVector::GetTotalMass () const
{
	return cZAIMass.GetMass(*this);
}
//____________________________________________________________________________
double IsotopicVector::GetMeanMolarMass () const
{
	return GetTotalMass() * 1e6 * AVOGADRO / GetActinidesComposition().GetSumOfAll();
}

//____________________________________________________________________________
std::vector< int > IsotopicVector::GetChemicalSpecies () const
{
	std::vector< int > tmp( fdata.size() );

	std::transform (
			fdata.begin() , fdata.end() , // input iterator
			tmp.begin()                 , // output iterator
			[]( const std::pair<ZAI,double> & zaiQ ) { return zaiQ.first.Z(); }
		);

	std::vector<int>::iterator end = std::unique ( tmp.begin() , tmp.end() );
	tmp.resize( end - tmp.begin() );

	return tmp;
}

//____________________________________________________________________________
double IsotopicVector::GetSumOfAll () const
{
	return std::accumulate(
			fdata.begin() , fdata.end() , // input iterator
			0.0 ,                         // first value for the sum
			[]( const double sum , const std::pair<ZAI,double> & zaiQ ) { return zaiQ.second + sum; }
		);
}

		
//____________________________________________________________________________
std::size_t IsotopicVector::size () const
	{ return fdata.size(); }
//____________________________________________________________________________
IsotopicVector::const_iterator IsotopicVector::begin () const
	{ return fdata.begin(); }
//____________________________________________________________________________
IsotopicVector::iterator IsotopicVector::begin ()
	{ return fdata.begin(); }
//____________________________________________________________________________
IsotopicVector::const_iterator IsotopicVector::end   () const
	{ return fdata.end(); }
//____________________________________________________________________________
IsotopicVector::iterator IsotopicVector::end   ()
	{ return fdata.end(); }
//____________________________________________________________________________
IsotopicVector::const_iterator IsotopicVector::find ( const ZAI & zai ) const
	{ return fdata.find(zai); }
//____________________________________________________________________________
IsotopicVector::iterator IsotopicVector::find ( const ZAI & zai )
	{ return fdata.find(zai); }


///// SETTER /////////////////////////////////////////////////////////////////
void IsotopicVector::Initialize ( double v )
{
	if ( v == 0 )
	{
		fdata.clear();
	}
	else
	{	
		iterator end = fdata.end();
		for ( iterator it=fdata.begin() ; it!=end ; ++it )
			{ it->second = v; }
	}
}
//____________________________________________________________________________
void IsotopicVector::Clear ()
	{ fdata.clear(); }
//____________________________________________________________________________
void IsotopicVector::Add ( const short int z , const short int a , const short int i , double q )
	{ fdata[ ZAI(z,a,i) ] += q; }
//____________________________________________________________________________
void IsotopicVector::Add ( const ZAI & zai , double q )
	{ fdata[zai] += q; }
//____________________________________________________________________________
void IsotopicVector::Add ( const std::pair<ZAI,double> & zaiQ )
	{ fdata.insert(zaiQ); }
//____________________________________________________________________________
void IsotopicVector::Add ( const IsotopicVector & a )
	{ fdata.insert( a.fdata.begin() , a.fdata.end() ); }
//____________________________________________________________________________
void IsotopicVector::Add ( const std::map<ZAI,double> & a )
	{ fdata.insert( a.begin() , a.end() ); }

//____________________________________________________________________________
void IsotopicVector::Remove ( const short int z , const short int a , const short int i , double q )
	{ Remove( ZAI(z,a,i) , q ); }
//____________________________________________________________________________
void IsotopicVector::Remove ( const ZAI & zai , double q )
{
	iterator it = fdata.find( zai );
	if ( it != fdata.end() )
	{
		it->second -= q;
		if ( it->second == 0 ) { fdata.erase(it); }
	}
}
//____________________________________________________________________________
void IsotopicVector::Remove ( const std::pair<ZAI,double> & zaiQ )
{
	iterator it = fdata.find(zaiQ.first);
	if ( it!= fdata.end() )
	{
		it->second -= zaiQ.second;
		if ( it->second == 0 ) { fdata.erase(it); }
	}
}
//____________________________________________________________________________
void IsotopicVector::Remove ( const IsotopicVector & a )
{
	const_iterator a_end = a.end();
	
	iterator f_end = fdata.end();
	iterator jt;
	for ( const_iterator it=a.begin() ; it!=a_end ; ++it )
	{
		jt = fdata.find( it->first );
		if ( jt != f_end )
		{
			jt->second -= it->second;
			if ( jt->second == 0 ) { fdata.erase(jt); }
		}
	}
}

//____________________________________________________________________________
void IsotopicVector::Multiply ( double f )
{
	iterator end = fdata.end();

	if ( f == 0 ) { fdata.clear(); }
	else
	{
		for ( iterator it=fdata.begin() ; it!=end ; ++it )
			{ it->second *= f; }
	}
}

//____________________________________________________________________________
void IsotopicVector::ApplyZAIThreshold ( int a )
{
	double quantity = 0.0;

	for ( iterator it = begin() ; it!=end() ; ++it )
	{
		if ( it->first.Z() < a )
		{
			quantity += it->second;
			it = fdata.erase(it);
		}
	}
	fdata.insert( std::pair<ZAI,double>( ZAI(0,0,0) , quantity ) );
}

// IO Methods
void IsotopicVector::Write ( std::string filename , cSecond time ) const
{
	std::ofstream f( filename.c_str() , std::ios::app );
	if ( !f.good() )
		{ std::cerr << "!!Warning!! !!!IsotopicVector!!! \n Can't open \"" << filename << "\"\n" << std::endl; }

	if ( time != -1 )
	{
		f << "Time" << time/cYear << std::endl;
	}
	std::transform(
			fdata.begin() , fdata.end() , // input iterator
			std::ostream_iterator<std::string>( f , "\n" ) , // output iterator
			[]( const std::pair<ZAI,double> & zaiQ ) {
				std::stringstream ss;
				ss << zaiQ.first.Z() << " " << zaiQ.first.A() << " " << zaiQ.first.I() << " ";
				ss << zaiQ.second;
				return ss.str();
			}
		);
}
//____________________________________________________________________________
void IsotopicVector::Print ( std::string o ) const
{
	std::cout << sPrint();
}
//____________________________________________________________________________
std::string IsotopicVector::sPrint() const
{
	std::stringstream ss;
	ss << "**************************\n";
	ss << "*Isotopic Vector Property*\n";
	ss << "**************************\n" << std::endl;

	ss << "*Isotopic Vector Quantity*" << std::endl;

	const_iterator end = fdata.end();
	for( const_iterator it = fdata.begin() ; it != end ; ++it )
	{
		ss << it->first.Z() << " ";
		ss << it->first.A() << " ";
		ss << it->first.I() << " ";
		ss << ": " << it->second;
		ss << "\n";
	}
	ss << endl;

	return ss.str();
}
//____________________________________________________________________________
void IsotopicVector::PrintList( std::string o ) const
{
	const_iterator end = fdata.end();
	for( const_iterator it = fdata.begin() ; it != end ; ++it )
	{
		std::cout << it->first.Z() << " ";
		std::cout << it->first.A() << " ";
		std::cout << it->first.I() << " ";
		std::cout << ": " << it->second;
		std::cout << "\n";
	}
	std::cout << endl;
}
//____________________________________________________________________________


///// EXTERN OPERATOR ////////////////////////////////////////////////////////
bool operator == ( const IsotopicVector & a , const IsotopicVector & b )
{
	bool r = true;
	
	IsotopicVector::const_iterator end = a.end();
	IsotopicVector::const_iterator it  = a.end();
	while ( it != end && r )
	{
		r = b[ it->first ] == it->second;
		++it; 
	}

	return r;
}
//____________________________________________________________________________
bool operator != ( const IsotopicVector & a , const IsotopicVector & b )
{
	bool r = true;
	
	IsotopicVector::const_iterator end = a.end();
	IsotopicVector::const_iterator it  = a.end();
	while ( it != end && r )
	{
		r = b[ it->first ] != it->second;
		++it; 
	}

	return r;
}
//____________________________________________________________________________
bool operator <  ( const IsotopicVector & a , const IsotopicVector & b )
{
	bool r = true;
	
	IsotopicVector::const_iterator end = a.end();
	IsotopicVector::const_iterator it  = a.end();
	while ( it != end && r )
	{
		r = b[ it->first ] <  it->second;
		++it; 
	}

	return r;
}
//____________________________________________________________________________
bool operator >  ( const IsotopicVector & a , const IsotopicVector & b )
{
	bool r = true;
	
	IsotopicVector::const_iterator end = a.end();
	IsotopicVector::const_iterator it  = a.end();
	while ( it != end && r )
	{
		r = b[ it->first ] >  it->second;
		++it; 
	}

	return r;
}
//____________________________________________________________________________
bool operator <= ( const IsotopicVector & a , const IsotopicVector & b )
{
	bool r = true;
	
	IsotopicVector::const_iterator end = a.end();
	IsotopicVector::const_iterator it  = a.end();
	while ( it != end && r )
	{
		r = b[ it->first ] <= it->second;
		++it; 
	}

	return r;
}
//____________________________________________________________________________
bool operator >= ( const IsotopicVector & a , const IsotopicVector & b )
{
	bool r = true;
	
	IsotopicVector::const_iterator end = a.end();
	IsotopicVector::const_iterator it  = a.end();
	while ( it != end && r )
	{
		r = b[ it->first ] >= it->second;
		++it; 
	}

	return r;
}
//____________________________________________________________________________
IsotopicVector operator / ( const IsotopicVector & a , double f )
{
	IsotopicVector tmp = a;
	IsotopicVector::iterator end = tmp.end();
	for ( IsotopicVector::iterator it=tmp.begin() ; it!=end ; ++it )
		{ it->second /= f; }

	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator / ( const ZAI & a , double f )
{
	IsotopicVector tmp;
	tmp.Add(a,1/f);
	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator * ( const IsotopicVector & a , double f )
{
	IsotopicVector tmp = a;
	IsotopicVector::iterator end = tmp.end();
	for ( IsotopicVector::iterator it=tmp.begin() ; it!=end ; ++it )
		{ it->second *= f; }

	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator * ( const ZAI & a , double f )
{
	IsotopicVector tmp;
	tmp.Add(a,f);
	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator * ( double f , const IsotopicVector & a )
{
	IsotopicVector tmp = a;
	IsotopicVector::iterator end = tmp.end();
	for ( IsotopicVector::iterator it=tmp.begin() ; it!=end ; ++it )
		{ it->second *= f; }

	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator * ( double f , const ZAI & a )
{
	IsotopicVector tmp;
	tmp.Add(a,f);
	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator + ( const IsotopicVector & a , const IsotopicVector & b )
{
	IsotopicVector tmp = a;
	tmp += b;
	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator - ( const IsotopicVector & a , const IsotopicVector & b )
{
	IsotopicVector tmp = a;
	tmp -= b;
	return tmp;
}
//____________________________________________________________________________
IsotopicVector operator * ( const IsotopicVector & a , const IsotopicVector & b )
{
	IsotopicVector tmp = a;
	tmp *= b;
	return tmp;
}

//____________________________________________________________________________
std::ostream & operator << ( std::ostream & os , const IsotopicVector & a )
{
	std::transform(
			a.begin() , a.end() ,                         // input iterator
			std::ostream_iterator<std::string>(os,"\n") , // output iterator
			[]( const std::pair<ZAI,double> & zaiQ )
				{
					std::stringstream ss; ss << zaiQ.first.Z() << " " << zaiQ.first.A() << " " << zaiQ.first.I();
					ss << " " << zaiQ.second;
					return ss.str();
				}
		);

	return os;
}

///// FUNCTIONS //////////////////////////////////////////////////////////////
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

///// NEEDED /////////////////////////////////////////////////////////////////
void IsotopicVector::ClearNeed() { fIsotopicQuantityNeeded.clear(); }
std::map<ZAI ,double> IsotopicVector::GetIsotopicQuantityNeeded () const
{
	return fIsotopicQuantityNeeded;
}
