//////////////////////////////////////////////////////////////////////////////
///// CLASS Array<T> /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///// CONSTRUCTOR ////////////////////////////////////////////////////////////
template <typename T>
Array<T>::Array () :
	fdata(0) , fxInit(0) , fxStep(0)
{ ; }

//____________________________________________________________________________
template <typename T>
Array<T>::Array ( const Array<T> & a ) :
	fdata(a.fdata) , fxInit(a.fxInit) , fxStep(a.fxStep)
{ ; }

//____________________________________________________________________________
template <typename T>
Array<T>::Array ( std::size_t size , T xInit , T xStep ) :
	fdata(size) , fxInit(xInit) , fxStep(xStep)
{ ; }

//____________________________________________________________________________
template <typename T>
Array<T>::Array ( const T * a , std::size_t size , T xInit , T xStep ) :
	fdata(a,size) , fxInit(xInit) , fxStep(xStep)
{ ; }

//____________________________________________________________________________
template <typename T>
Array<T>::Array ( const std::valarray<T> & a , T xInit , T xStep ) :
	fdata(a) , fxInit(xInit) , fxStep(xStep)
{ ; }

//____________________________________________________________________________
template <typename T>
Array<T>::Array ( const std::vector<T> & a , T xInit , T xStep ) :
	fdata(&a[0],a.size()) , fxInit(xInit) , fxStep(xStep)
{ ; }

//____________________________________________________________________________
template <typename T> template <typename F>
Array<T>::Array ( std::size_t size , F f , T xInit , T xStep ) :
	fdata(size) , fxInit(xInit) , fxStep(xStep)
{
	std::generate( begin() , end() , f );
}


///// DESTRUCTOR /////////////////////////////////////////////////////////////
template <typename T>
Array<T>::~Array ()
{ ; }

///// OPERATOR ///////////////////////////////////////////////////////////////
template <typename T>
Array<T> & Array<T>::operator = ( const Array<T> & a )
{
	set( a );
	return *this;
}

//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator += ( const Array<T> & a )
{
	//if ( fdata.size() != a.fdata.size() )
	//	{ throw std::length_error(__PRETTY_FUNCTION__ + ": two arrays have not the same size"); }

	fdata += a.fdata;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator -= ( const Array<T> & a )
{
	//if ( fdata.size() != a.fdata.size() )
	//	{ throw std::length_error(__PRETTY_FUNCTION__ + ": two arrays have not the same size"); }
	
	fdata -= a.fdata;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator *= ( const Array<T> & a )
{
	//if ( fdata.size() != a.fdata.size() )
	//	{ throw std::length_error(__PRETTY_FUNCTION__ + ": two arrays have not the same size"); }
	
	fdata *= a.fdata;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator /= ( const Array<T> & a )
{
	//if ( fdata.size() != a.fdata.size() )
	//	{ throw std::length_error(__PRETTY_FUNCTION__ + ": two arrays have not the same size"); }
	
	fdata /= a.fdata;
	return *this;
}

//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator += ( const T & v )
{
	fdata += v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator -= ( const T & v )
{
	fdata -= v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator *= ( const T & v )
{
	fdata *= v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array<T> & Array<T>::operator /= ( const T & v )
{
	//if ( v == 0 )
	//	{ throw std::invalid_argument(__PRETTY_FUNCTION__ + ": you try something wrong (you can't divide by zero !!!)"); }

	fdata /= v;
	return *this;
}

//____________________________________________________________________________
template <typename T>
Array<T> Array<T>::operator + () const
{
	return Array<T>( +fdata , fxInit , fxStep );
}
//____________________________________________________________________________
template <typename T>
Array<T> Array<T>::operator - () const
{
	return Array<T>( -fdata , fxInit , fxStep );
}

//____________________________________________________________________________
template <typename T>
T & Array<T>::operator [] ( std::size_t i )
{
	//if ( i < 0 || i > fdata.size() )
	//	{ throw std::domain_error(__PRETTY_FUNCTION__ + ": out of bounds"); }

	return fdata[i];
}
//____________________________________________________________________________
template <typename T>
T Array<T>::operator [] ( std::size_t i ) const
{
	//if ( i < 0 || i > fdata.size() )
	//	{ throw std::domain_error(__PRETTY_FUNCTION__ + ": out of bounds"); }

	return fdata[i];
}


///// GETTER /////////////////////////////////////////////////////////////////
template <typename T>
std::size_t Array<T>::size () const
{ return fdata.size(); }

//____________________________________________________________________________
template <typename T>
T Array<T>::at ( std::size_t i ) const
{
	//if ( i < 0 || i > fdata.size() )
	//	{ throw std::domain_error(__PRETTY_FUNCTION__ + ": out of bounds"); }
	
	return fdata[i];
}

//____________________________________________________________________________
template <typename T>
std::size_t Array<T>::getBin ( T x ) const
{
	return (std::size_t)( (x - fxInit)/fxStep );
}
//____________________________________________________________________________
template <typename T>
T Array<T>::getX ( std::size_t i ) const
{
	return fxInit + fxStep*(T)i;
}

//____________________________________________________________________________
template <typename T>
T Array<T>::xBegin () const
{
	return fxInit;
}
//____________________________________________________________________________
template <typename T>
T Array<T>::xEnd () const
{
	return fxInit + fxStep*fdata.size();
}
//____________________________________________________________________________
template <typename T>
T Array<T>::xStep () const
{
	return fxStep;
}

//____________________________________________________________________________
template <typename T>
T * Array<T>::begin () { return &fdata[0]; }
//____________________________________________________________________________
template <typename T>
T * Array<T>::end   () { return &fdata[ fdata.size()-1 ]+1; }

///// METHOD /////////////////////////////////////////////////////////////////
template <typename T>
T Array<T>::min () const
{ return fdata.min(); }
//____________________________________________________________________________
template <typename T>
T Array<T>::max () const
{ return fdata.max(); }
//____________________________________________________________________________
template <typename T>
T Array<T>::sum () const
{ return fdata.sum(); }

//____________________________________________________________________________
template <typename T>
std::valarray<T> & Array<T>::data ()
{ return fdata; }
//____________________________________________________________________________
template <typename T>
std::valarray<T> Array<T>::data () const
{ return fdata; }

//____________________________________________________________________________
template <typename T> template <typename F>
Array<T> Array<T>::apply ( F f ) const
{ return Array<T>( fdata.apply(f) , fxInit , fxStep ); }

//____________________________________________________________________________
template <typename T>
T Array<T>::eval ( T x ) const
{
	/* interpolation linéaire :
		$$y = y_a + \frac{y_b - y_a}{x_b - x_a}(x - x_a)$$
	*/
	std::size_t i = getBin(x); // indice de x_a
	const T xa = getX(i) , xb = getX(i+1);

	T y = fdata[i] + ( fdata[i+1] - fdata[i] )/( xb - xa ) * (x - xa);
	return y;
}

//____________________________________________________________________________
template <typename T>
void Array<T>::set ( const Array<T> & a )
{
	if ( this != &a )
	{
		fdata = a.fdata;
		fxInit = a.fxInit;
		fxStep = a.fxStep;
	}
}


///// OPERATOR ///////////////////////////////////////////////////////////////
template <typename T>
Array<bool> operator < ( const Array<T> & a , const Array<T> & b )
{
	return Array<bool> ( a.data() < b.data() , a.xBegin() < b.xBegin() , a.xStep() < b.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<bool> operator >  ( const Array<T> & a , const Array<T> & b )
{
	return Array<bool> ( a.data() > b.data() , a.xBegin() > b.xBegin() , a.xStep() > b.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<bool> operator == ( const Array<T> & a , const Array<T> & b )
{
	return Array<bool> ( a.data() == b.data() , a.xBegin() == b.xBegin() , a.xStep() == b.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<bool> operator <=  ( const Array<T> & a , const Array<T> & b )
{
	return Array<bool> ( a.data() <= b.data() , a.xBegin() <= b.xBegin() , a.xStep() <= b.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<bool> operator >=  ( const Array<T> & a , const Array<T> & b )
{
	return Array<bool> ( a.data() >= b.data() , a.xBegin() >= b.xBegin() , a.xStep() >= b.xStep() );
}


//____________________________________________________________________________
template <typename T>
Array<T> operator + ( const Array<T> & a , const Array<T> & b )
{
	return Array<T>( a.data()+b.data() , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator - ( const Array<T> & a , const Array<T> & b )
{
	return Array<T>( a.data()-b.data() , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator * ( const Array<T> & a , const Array<T> & b )
{
	return Array<T>( a.data()*b.data() , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator / ( const Array<T> & a , const Array<T> & b )
{
	return Array<T>( a.data()/b.data() , a.xBegin() , a.xStep() );
}

//____________________________________________________________________________
template <typename T>
Array<T> operator + ( const T & v , const Array<T> & a )
{
	return Array<T>( v+a.data() , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator - ( const T & v , const Array<T> & a )
{
	return Array<T>( v-a.data() , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator * ( const T & v , const Array<T> & a )
{
	return Array<T>( v*a.data() , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator / ( const T & v , const Array<T> & a )
{
	return Array<T>( v/a.data() , a.xBegin() , a.xStep() );
}

//____________________________________________________________________________
template <typename T>
Array<T> operator + ( const Array<T> & a , const T & v )
{
	return Array<T>( a.data()+v , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator - ( const Array<T> & a , const T & v )
{
	return Array<T>( a.data()-v , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator * ( const Array<T> & a , const T & v )
{
	return Array<T>( a.data()*v , a.xBegin() , a.xStep() );
}
//____________________________________________________________________________
template <typename T>
Array<T> operator / ( const Array<T> & a , const T & v )
{
	return Array<T>( a.data()/v , a.xBegin() , a.xStep() );
}
