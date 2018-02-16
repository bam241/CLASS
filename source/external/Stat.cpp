//////////////////////////////////////////////////////////////////////////////
///// CLASS Stat<T> //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///// CONSTRUCTOR ////////////////////////////////////////////////////////////
template <typename T>
Stat<T>::Stat () :
	fn(0) , fmin(0) , fmax(0) , fsum(0) , fsum2(0)
{ ; }

//____________________________________________________________________________
template <typename T>
Stat<T>::Stat ( const Stat<T> & a ) :
	fn(a.fn) , min(a.fmin) , fmax(a.fmax) , fsum(a.fsum) , fsum2(a.fsum2)
{ ; }

//____________________________________________________________________________
template <typename T> template <typename ITR>
Stat<T>::Stat ( const ITR first , const ITR last ) :
	fn(0) , fmin(*first) , fmax(*first) , fsum(0) , fsum2(0)
{
	for ( ITR it = first ; it != last ; ++it ) // do all in one loop
	{
		if ( *it > fmax ) { fmax = *it; }
		if ( *it < fmin ) { fmin = *it; }

		fsum  += *it;
		fsum2 += (*it) * (*it);
		++fn;
	}
}

///// DESTRUCTOR /////////////////////////////////////////////////////////////
template <typename T>
Stat<T>::~Stat ()
{ ; }

///// OPERATOR ///////////////////////////////////////////////////////////////
template <typename T>
Stat<T> & Stat<T>::operator = ( const Stat<T> & a )
{
	if ( this != &a )
	{
		fn = a.fn;
		fmin = a.fmin; fmax  = a.fmax;
		fsum = a.fsum; fsum2 = a.fsum2;
	}
	return *this;
}

//____________________________________________________________________________
template <typename T>
Stat<T> & Stat<T>::operator += ( const T t )
{ this->add( t ); }

///// GETTER /////////////////////////////////////////////////////////////////
template <typename T> std::size_t Stat<T>::n () const { return fn; }

//____________________________________________________________________________
template <typename T> T Stat<T>::min () const  { return fmin; }

//____________________________________________________________________________
template <typename T> T Stat<T>::max () const  { return fmax; }

//____________________________________________________________________________
template <typename T> T Stat<T>::amp () const  { return fmax - fmin; }

//____________________________________________________________________________
template <typename T> T Stat<T>::avg () const  { return fsum/fn;     }

//____________________________________________________________________________
template <typename T> T Stat<T>::rms () const  { return std::sqrt( (fsum*fsum)/fn - fsum2/fn ); }

///// SETTER /////////////////////////////////////////////////////////////////
template <typename T>
void Stat<T>::add ( T val )
{
	fn   += 1;
	fmin = std::min( fmin , val ); fsum  += val;
	fmax = std::max( fmax , val ); fsum2 += val*val;
}
