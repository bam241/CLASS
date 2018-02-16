//////////////////////////////////////////////////////////////////////////////
///// CLASS Array<T> /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///// CONSTRUCTOR ////////////////////////////////////////////////////////////
template <typename T>
Array<T>::Array () :
	fdata() , fxInit(0) , fxStep(0)
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
	std::generate( &fdata[0] , &fdata[size] , f );
}


///// DESTRUCTOR /////////////////////////////////////////////////////////////
template <typename T>
Array<T>::~Array ()
{ ; }

///// OPERATOR ///////////////////////////////////////////////////////////////
template <typename T>
Array & Array<T>::operator = ( const Array<T> & a )
{
	set( a );
	return *this;
}

//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator += ( const Array<T> & a )
{
	fdata += a.fdata;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator -= ( const Array<T> & a )
{
	fdata -= a.fdata;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator *= ( const Array<T> & a )
{
	fdata *= a.fdata;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator /= ( const Array<T> & a )
{
	fdata /= a.fdata;
	return *this;
}

//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator += ( const T & v )
{
	fdata += v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator -= ( const T & v )
{
	fdata -= v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator *= ( const T & v )
{
	fdata *= v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Array & Array<T>::operator /= ( const T & v )
{
	fdata /= v;
	return *this;
}

//____________________________________________________________________________
template <typename T>
Array Array<T>::operator + () const
{
	return Array<T>( +fdata , fxInit , fxStep );
}
//____________________________________________________________________________
template <typename T>
Array Array<T>::operator - () const
{
	return Array<T>( -fdata , fxInit , fxStep );
}

//____________________________________________________________________________
template <typename T>
T & Array<T>::operator [] ( std::size_t i )
{
	return fdata[i];
}
//____________________________________________________________________________
template <typename T>
T Array<T>::operator [] ( std::size_t i ) const
{
	return fdata[i];
}


///// GETTER /////////////////////////////////////////////////////////////////
template <typename T>
std::size_t Array<T>::size () const
{ return fdata.size(); }

//____________________________________________________________________________
template <typename T>
T Array<T>::at ( std::size_t i ) const
{ return fdata[i]; }

//____________________________________________________________________________
template <typename T>
std::size_t Array<T>::getBin ( T x ) const
{
	return (std::size_t)(x - fxInit)/fxStep;
}
//____________________________________________________________________________
template <typename T>
std::size_t Array<T>::getX ( std::size_t i ) const
{
	return fxInit + (i+0.0)*fxStep;
}

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
