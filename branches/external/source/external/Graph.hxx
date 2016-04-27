#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <valarray>
#include <algorithm>

#include "Array.hxx"

template <typename T>
class Graph : public Array<T>
{
	public:
		// CONSTRUCTOR
		Graph ();
		Graph ( const Graph<T> & );

		Graph ( std::size_t );
		
		Graph ( std::size_t , const T * );
		Graph ( std::size_t , const std::valarray<T> & );
		Graph ( std::size_t , const std::vector<T>   & );
		Graph ( const T * , std::size_t  , const T * );
		Graph ( const T * , std::size_t  , const std::valarray<T> & );
		Graph ( const T * , std::size_t  , const std::vector<T>   & );
		Graph ( const std::valarray<T> & , const T * );
		Graph ( const std::valarray<T> & , const std::valarray<T> & );
		Graph ( const std::valarray<T> & , const std::vector<T>   & );
		Graph ( const std::vector<T>   & , const T * );
		Graph ( const std::vector<T>   & , const std::valarray<T> & );
		Graph ( const std::vector<T>   & , const std::vector<T>   & );
		template <typename F>
		Graph ( std::size_t , F , const T * );
		template <typename F>
		Graph ( std::size_t , F , const std::valarray<T> & );
		template <typename F>
		Graph ( std::size_t , F , const std::vector<T>   & );

		// DESTRUCTOR
		virtual ~Graph ();

		// OPERATOR
		virtual Graph<T> & operator = ( const Graph<T> & );

		Graph & operator += ( const Array<T> & );
		Graph & operator -= ( const Array<T> & );
		Graph & operator *= ( const Array<T> & );
		Graph & operator /= ( const Array<T> & );

		Graph & operator += ( const T & );
		Graph & operator -= ( const T & );
		Graph & operator *= ( const T & );
		Graph & operator /= ( const T & );

		Graph   operator + () const;
		Graph   operator - () const;

		// GETTER
		virtual std::size_t getBin ( T )           const;
		virtual T           getX   ( std::size_t ) const;

		virtual T xBegin () const;
		virtual T xEnd   () const;

		// METHOD
		std::valarray<T> & time ();
		std::valarray<T>   time () const;

		template <typename F>
		Graph<T> apply ( F ) const;

		T eval ( T ) const;

	protected:
		virtual void set ( const Graph<T> & );

	private:
		std::valarray<T> ftime;

};

template <typename T>
Graph<bool> operator <  ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<bool> operator >  ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<bool> operator == ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<bool> operator <= ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<bool> operator >= ( const Graph<T> & , const Graph<T> & );

template <typename T>
Graph<T> operator + ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<T> operator - ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<T> operator * ( const Graph<T> & , const Graph<T> & );
template <typename T>
Graph<T> operator / ( const Graph<T> & , const Graph<T> & );

template <typename T>
Graph<T> operator + ( const T & , const Graph<T> & );
template <typename T>
Graph<T> operator - ( const T & , const Graph<T> & );
template <typename T>
Graph<T> operator * ( const T & , const Graph<T> & );
template <typename T>
Graph<T> operator / ( const T & , const Graph<T> & );

template <typename T>
Graph<T> operator + ( const Graph<T> & , const T & );
template <typename T>
Graph<T> operator - ( const Graph<T> & , const T & );
template <typename T>
Graph<T> operator * ( const Graph<T> & , const T & );
template <typename T>
Graph<T> operator / ( const Graph<T> & , const T & );

#include "Graph.cpp"

#endif
