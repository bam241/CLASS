//////////////////////////////////////////////////////////////////////////////
///// CLASS Graph<T> /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

///// CONSTRUCTOR ////////////////////////////////////////////////////////////
template <typename T>
Graph<T>::Graph () :
	Array<T>() , ftime()
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const Graph<T> & a ) :
	Array<T>(a) , ftime(a.ftime)
{ ; }

//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( std::size_t size ) :
	Array<T>(size) , ftime(size)
{ ; }

//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( std::size_t size , const T * time ) :
	Array<T>(size) , ftime(time,size)
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( std::size_t size , const std::valarray<T> & time ) :
	Array<T>(size) , ftime(&time[0],size)
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( std::size_t size , const std::vector<T> & time ) :
	Array<T>(size) , ftime(&time[0],size)
{ ; }

//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const T * a , std::size_t size , const T * time ) :
	Array<T>(a,size) , ftime(time,size)
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const T * a , std::size_t size , const std::valarray<T> & time ) :
	Array<T>(a,size) , ftime(&time[0],size)
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const T * a , std::size_t size , const std::vector<T> & time ) :
	Array<T>(a,size) , ftime(&time[0],size)
{ ; }

//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::valarray<T> & a , const T * time ) :
	Array<T>(a) , ftime(time,a.size())
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::valarray<T> & a , const std::valarray<T> & time ) :
	Array<T>(a) , ftime(&time[0],a.size())
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::valarray<T> & a , const std::vector<T> & time ) :
	Array<T>(a) , ftime(&time[0],a.size())
{ ; }

//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::vector<T> & a , const T * time ) :
	Array<T>(a) , ftime(time,a.size())
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::vector<T> & a , const std::valarray<T> & time ) :
	Array<T>(a) , ftime(&time[0],a.size())
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::vector<T> & a , const std::vector<T> & time ) :
	Array<T>(a) , ftime(&time[0],a.size())
{ ; }

//____________________________________________________________________________
template <typename T> template <typename F>
Graph<T>::Graph ( std::size_t size , F f , const T * time ) :
	Array<T>(size,f) , ftime(time,size)
{ ; }
//____________________________________________________________________________
template <typename T> template <typename F>
Graph<T>::Graph ( std::size_t size , F f , const std::valarray<T> & time ) :
	Array<T>(size,f) , ftime(&time[0],size)
{ ; }
//____________________________________________________________________________
template <typename T> template <typename F>
Graph<T>::Graph ( std::size_t size , F f , const std::vector<T> & time ) :
	Array<T>(size,f) , ftime(&time[0],size)
{ ; }

///// DESTRUCTOR /////////////////////////////////////////////////////////////
template <typename T>
Graph<T>::~Graph ()
{ ; }

///// OPERATOR ///////////////////////////////////////////////////////////////
template <typename T>
Graph<T> & Graph<T>::operator = ( const Graph<T> & a )
{
	set( a );
	return *this;
}

//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator += ( const Array<T> & a )
{
	this->data() += a.data();
	return *this;
}
//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator -= ( const Array<T> & a )
{
	this->data() -= a.data();
	return *this;
}
//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator *= ( const Array<T> & a )
{
	this->data() *= a.data();
	return *this;
}
//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator /= ( const Array<T> & a )
{
	this->data() /= a.data();
	return *this;
}

//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator += ( const T & v )
{
	this->data() += v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator -= ( const T & v )
{
	this->data() -= v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator *= ( const T & v )
{
	this->data() *= v;
	return *this;
}
//____________________________________________________________________________
template <typename T>
Graph<T> & Graph<T>::operator /= ( const T & v )
{
	this->data() /= v;
	return *this;
}

//____________________________________________________________________________
template <typename T>
Graph<T> Graph<T>::operator + () const
{
	return Graph<T>( +(this->data()) , &ftime[0] );
}
//____________________________________________________________________________
template <typename T>
Graph<T> Graph<T>::operator - () const
{
	return Graph<T>( -(this->data()) , &ftime[0] );
}

///// GETTER /////////////////////////////////////////////////////////////////
template <typename T>
std::size_t Graph<T>::getBin ( T x ) const
{
	// /!\ WARNING : be carefull  with x
	std::size_t start=0 , end=this->size() , mid=(start+end)/2;
	bool find = false;
	do
	{
		mid=(start+end)/2;
		if ( this->at(mid) == x )
			{ find = true; }
		else
		{
			if ( x > this->at(mid) )
				{ start = mid+1; }
			else
				{ end = mid -1; }
		}
	} while ( !find && start < end );

	return mid;
}
//____________________________________________________________________________
template <typename T>
T Graph<T>::getX ( std::size_t i ) const
{
	return ftime[i];
}

//____________________________________________________________________________
template <typename T>
T Graph<T>::xBegin () const
{ return ftime[0]; }
//____________________________________________________________________________
template <typename T>
T Graph<T>::xEnd () const
{ return ftime[ftime.size()]; }

///// METHOD /////////////////////////////////////////////////////////////////
template <typename T>
std::valarray<T> & Graph<T>::time ()
{ return ftime; }
//____________________________________________________________________________
template <typename T>
std::valarray<T> Graph<T>::time () const
{ return ftime; }

//____________________________________________________________________________
template <typename T> template <typename F>
Graph<T> Graph<T>::apply ( F f ) const
{ return Graph<T>( this->data().apply(f) , &ftime[0] ); }

//____________________________________________________________________________
template <typename T>
T Graph<T>::eval ( T x ) const
{
	std::size_t i = getBin ( x );

	const T xa = getX(i) , xb = getX(i+1);

	T y = this->at(i) + ( this->at(i+1) - this->at(i) )/(xb - xa) * (x - xa);

	return y;
}

//____________________________________________________________________________
template <typename T>
void Graph<T>::set ( const Graph<T> &  a )
{
	if ( this != &a )
	{
		Array<T>::set ( a );
		ftime = a.ftime;
	}
}


///// OPERATOR ///////////////////////////////////////////////////////////////
template <typename T>
Graph<bool> operator < ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<bool> ( a.data() < b.data() , a.time() < b.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<bool> operator > ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<bool> ( a.data() > b.data() , a.time() > b.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<bool> operator == ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<bool> ( a.data() == b.data() , a.time() == b.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<bool> operator <= ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<bool> ( a.data() <= b.data() , a.time() <= b.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<bool> operator >= ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<bool> ( a.data() >= b.data() , a.time() >= b.time() );
}

//____________________________________________________________________________
template <typename T>
Graph<T> operator + ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<T> ( a.data() + b.data() , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator - ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<T> ( a.data() - b.data() , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator * ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<T> ( a.data() * b.data() , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator / ( const Graph<T> & a , const Graph<T> & b )
{
	return Graph<T> ( a.data() / b.data() , a.time() );
}

//____________________________________________________________________________
template <typename T>
Graph<T> operator + ( const T & v , const Graph<T> & a )
{
	return Graph<T> ( v + a.data() , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator - ( const T & v , const Graph<T> & a )
{
	return Graph<T> ( v - a.data() , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator * ( const T & v , const Graph<T> & a )
{
	return Graph<T> ( v * a.data() , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator / ( const T & v , const Graph<T> & a )
{
	return Graph<T> ( v / a.data() , a.time() );
}

//____________________________________________________________________________
template <typename T>
Graph<T> operator + ( const Graph<T> & a , const T & v )
{
	return Graph<T>( a.data()+v , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator - ( const Graph<T> & a , const T & v )
{
	return Graph<T>( a.data()-v , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator * ( const Graph<T> & a , const T & v )
{
	return Graph<T>( a.data()*v , a.time() );
}
//____________________________________________________________________________
template <typename T>
Graph<T> operator / ( const Graph<T> & a , const T & v )
{
	return Graph<T>( a.data()/v , a.time() );
}
