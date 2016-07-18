#include "Graph.hxx"

ClassImp(Graph)

///// CONSTRUCTOR ////////////////////////////////////////////////////////////
Graph::Graph ()
{ ; }
//____________________________________________________________________________
Graph::Graph ( const Graph & a ) :
	fx(a.fx) , fy(a.fy)
{ ; }
//____________________________________________________________________________
Graph::Graph ( std::size_t n ) :
	fx(n) , fy(n)
{ ; }
//____________________________________________________________________________
Graph::Graph ( std::size_t n , const double * x , const double * y ) :
	fx(x,n) , fy(n)
{
	if ( y != nullptr )
		{ std::copy( y , y+n , std::begin(fy) ); }
}

//____________________________________________________________________________
Graph::Graph ( const TGraph & a ) :
	fx( a.GetX() , a.GetN() ) , fy( a.GetY() , a.GetN() )
{ ; }
//____________________________________________________________________________
Graph::Graph ( const TGraph * a ) :
	fx( a->GetX() , a->GetN() ) , fy( a->GetY() , a->GetN() )
{ ; }

///// DESTRUCTOR /////////////////////////////////////////////////////////////
Graph::~Graph ()
{ ; }

///// OPERATOR ///////////////////////////////////////////////////////////////
Graph & Graph::operator = ( const Graph & a )
{
	if ( this != &a )
	{
		fx = a.fx;
		fy = a.fy;
	}
	return *this;
}

//____________________________________________________________________________
void Graph::operator = ( const TGraph & a )
{
	fx.resize( a.GetN() ); fy.resize( a.GetN() );
	std::copy( a.GetX() , a.GetX() + a.GetN() , std::begin(fx) );
	std::copy( a.GetY() , a.GetY() + a.GetN() , std::begin(fy) );
}
//____________________________________________________________________________
void Graph::operator = ( const TGraph * a )
{
	fx.resize( a->GetN() ); fy.resize( a->GetN() );
	std::copy( a->GetX() , a->GetX() + a->GetN() , std::begin(fx) );
	std::copy( a->GetY() , a->GetY() + a->GetN() , std::begin(fy) );
}

//____________________________________________________________________________
Graph::operator TGraph () const
{
	return TGraph( fx.size() , &fx[0] , &fy[0] );
}

//____________________________________________________________________________
Graph & Graph::operator += ( const Graph & a )
{
	fy += a.fy;
	return *this;
}
//____________________________________________________________________________
Graph & Graph::operator -= ( const Graph & a )
{
	fy -= a.fy;
	return *this;
}
//____________________________________________________________________________
Graph & Graph::operator *= ( const Graph & a )
{
	fy *= a.fy;
	return *this;
}

//____________________________________________________________________________
const double & Graph::operator [] ( std::size_t i ) const { return fy[i]; }
//____________________________________________________________________________
double &       Graph::operator [] ( std::size_t i )       { return fy[i]; }

///// GETTER /////////////////////////////////////////////////////////////////
std::size_t    Graph::size   () const { return fx.size();      }
//____________________________________________________________________________
const double * Graph::Xbegin () const { return std::begin(fx); }
//____________________________________________________________________________
double *       Graph::Xbegin ()       { return std::begin(fx); }
//____________________________________________________________________________
const double * Graph::Ybegin () const { return std::begin(fy); }
//____________________________________________________________________________
double *       Graph::Ybegin ()       { return std::begin(fy); }
//____________________________________________________________________________
const double * Graph::Xend   () const { return std::end(fx);   }
//____________________________________________________________________________
double *       Graph::Xend   ()       { return std::end(fx);   }
//____________________________________________________________________________
const double * Graph::Yend   () const { return std::end(fy);   }
//____________________________________________________________________________
double *       Graph::Yend   ()       { return std::end(fy);   }

//____________________________________________________________________________
const double * Graph::GetX () const { return std::begin(fx); }
//____________________________________________________________________________
double * Graph::GetX () { return std::begin(fx); }
//____________________________________________________________________________
const double * Graph::GetY () const { return std::begin(fy); }
//____________________________________________________________________________
double * Graph::GetY () { return std::begin(fy); }

//____________________________________________________________________________
double Graph::GetX ( std::size_t i ) const
{ fy[i]; }
//____________________________________________________________________________
double Graph::GetY ( std::size_t i ) const
{ fy[i]; }
//____________________________________________________________________________
double Graph::GetY ( double x ) const
{ Eval(x); }
//____________________________________________________________________________
std::size_t Graph::GetN () const
{ return fx.size(); }
//____________________________________________________________________________
int Graph::GetPoint ( int i , double & x , double & y ) const
{
	if ( i<0 || i>fx.size() ) { return -1; }

	x = fx[i];
	y = fy[i];

	return i;
}

///// SETTER /////////////////////////////////////////////////////////////////
void Graph::resize ( std::size_t size , double x , double y )
{ fx.resize(size,x); fy.resize(size,y); }
//____________________________________________________________________________
void Graph::SetX ( std::size_t i , double x )
{ fx[i] = x; }
//____________________________________________________________________________
void Graph::SetY ( std::size_t i , double y )
{ fy[i] = y; }
//____________________________________________________________________________
void Graph::SetPoint ( std::size_t i , double x , double y )
{
	fx[i] = x; fy[i] = y;
}

///// METHOD /////////////////////////////////////////////////////////////////
double Graph::Eval ( double x ) const
{
	if ( x < fx[0] )             { return fy[0]; }             // if x is lower than min value
	if ( x > fx[ fx.size()-1 ] ) { return fy[ fy.size()-1 ]; } // if x is greater than max value

	const double * it = std::lower_bound( std::begin(fx) , std::end(fx) , x );

	if ( *it == x ) { return fy[ it - std::begin(fx) ]; } // if x is an actual value

	const std::size_t i = it - std::begin(fx); // interpolation between [i-1;i]
	// $y = y_{i-1} + ( x - x_{i-1}) \frac{ y_i - y_{i-1} }{ x_i - x_{i-1} }$
	return fy[i-1] + ( x - fx[i-1] ) * ( fy[i] - fy[i-1] )/( fx[i] - fx[i-1] );
}
//____________________________________________________________________________
double Graph::ExpEval ( double x ) const
{
	if ( x < fx[0] )             { return fy[0]; }             // if x is lower than min value
	if ( x > fx[ fx.size()-1 ] ) { return fy[ fy.size()-1 ]; } // if x is greater than max value

	const double * it = std::lower_bound( std::begin(fx) , std::end(fx) , x );

	if ( *it == x ) { return fy[ it - std::begin(fx) ]; } // if x is an actual value

	const std::size_t i = it - std::begin(fx); // exponential interpolation between [i-1;i]
	// $y = A e^{\lambda t$ with $\lambda = \frac{ln(y_1) - ln(y_2)}{t_2 - t_1}$ and $A=y_2 e^{\lambda t_2}$
	double lambda = std::log( fy[i-1]/fy[i] )/( fx[i] - fx[i-1] );
	double A = fy[i] * std::exp( lambda * fx[i] ); 
	return A*std::exp( lambda*x );
}

//____________________________________________________________________________
void Graph::Sort ()
{
	// besoin de directement faire un algo de tri sur une structure de tableau, plutôt que créer un tableau de structure pour le trier puis le transformer
	std::vector< std::pair<double,double> > tmp ( fx.size() ); // array with pair of <x,y>
	std::transform(
			std::begin(fx) , std::end(fx) , std::begin(fy) , // input iterators
			tmp.begin() , // output iterator
			[]( double x , double y ) { return std::pair<double,double>(x,y); }
		);
	
	std::sort(
			tmp.begin() , tmp.end() , // input iterator
			[]( const std::pair<double,double> & a , const std::pair<double,double> & b ) // comparator on x
				{ return a.first < b.first; }
		);

	std::transform(
			tmp.begin() , tmp.end() , // input iterator
			std::begin(fx) , // output iterator
			[]( std::pair<double,double> & p ) { return p.first; }
		);
	std::transform(
			tmp.begin() , tmp.end() , // input iterator
			std::begin(fy) , // output iterator
			[]( std::pair<double,double> & p ) { return p.second; }
		);
}

///// EXTERNAL OPERATOR //////////////////////////////////////////////////////
Graph operator + ( const Graph & a , const Graph & b )
{
	Graph tmp (a);
	for ( std::size_t i=0 ; i<tmp.size() ; ++i )
		{ tmp[i] += b[i]; }
	return tmp;
}
//____________________________________________________________________________
Graph operator - ( const Graph & a , const Graph & b )
{
	Graph tmp (a);
	for ( std::size_t i=0 ; i<tmp.size() ; ++i )
		{ tmp[i] -= b[i]; }
	return tmp;
}
//____________________________________________________________________________
Graph operator * ( const Graph & a , const Graph & b )
{
	Graph tmp (a);
	for ( std::size_t i=0 ; i<tmp.size() ; ++i )
		{ tmp[i] *= b[i]; }
	return tmp;
}
//____________________________________________________________________________
Graph operator * ( double x , const Graph & a )
{
	Graph tmp (a);
	for ( std::size_t i=0 ; i<tmp.size() ; ++i )
		{ tmp[i] *= x; }
	return tmp;
}
//____________________________________________________________________________
Graph operator * ( const Graph & a , double x )
{
	Graph tmp (a);
	for ( std::size_t i=0 ; i<tmp.size() ; ++i )
		{ tmp[i] *= x; }
	return tmp;
}
//____________________________________________________________________________
Graph operator / ( const Graph & a , double x )
{
	Graph tmp (a);
	for ( std::size_t i=0 ; i<tmp.size() ; ++i )
		{ tmp[i] /= x; }
	return tmp;
}
