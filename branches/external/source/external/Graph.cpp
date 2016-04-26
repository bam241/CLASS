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
Graph<T>::Graph ( const T * a , std::size_t size , const T * time ) :
	Array<T>(a,size) , ftime(time,size)
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::valarray<T> & a , const T * time ) :
	Array<T>(a) , ftime(time,a.size())
{ ; }
//____________________________________________________________________________
template <typename T>
Graph<T>::Graph ( const std::vector<T> & a , const T * time ) :
	Array<T>(a) , ftime(time,a.size())
{ ; }

//____________________________________________________________________________
template <typename T> template <typename F>
Graph<T>::Graph ( std::size_t size , F f , const T * time ) :
	Array<T>(size,f) , ftime(time,size)
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
std::size_t Graph<T>::getBin ( T x ) const
{
	std::size_t start=0 , end=size() , mid=(start+end)/2;
	bool find = false;
	do
	{
		mid=(start+end)/2;
		if ( at(mid) == x )
			{ find = true; }
		else
		{
			if ( x > at(mid) )
				{ start = mid+1; }
			else
				{ end = mid -1; }
		}

	} while ( !find && start <= end )

	return mid;
}
//____________________________________________________________________________
template <typename T>
T Graph<T>::getX ( std::size_t i ) const
{
	return ftime[i];
}


///// METHOD /////////////////////////////////////////////////////////////////
template <typename T>
T Graph<T>::eval ( T x ) const
{
	// i = find in sorting array
	std::size_t i = getBin ( x );


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
