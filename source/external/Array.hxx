#ifndef _ARRAY_H
#define _ARRAY_H

#include <iostream>
#include <valarray>
#include <algorithm>

template<typename T>
class Array
{
	public:
		// CONSTRUCTOR
		Array ();
		Array ( const Array<T> & );
		Array ( std::size_t , T xInit=0 , T xStep=1 );
		Array ( const T * , std::size_t  , T xInit=0 , T xStep=1 );
		Array ( const std::valarray<T> & , T xInit=0 , T xStep=1 );
		Array ( const std::vector<T>   & , T xInit=0 , T xStep=1 );
		template <typename F>
		Array ( std::size_t , F , T xInit=0 , T xStep=1 );

		// DESTRUCTOR
		virtual ~Array ();

		// OPERATOR
		virtual Array & operator = ( const Array<T> & );

		Array & operator += ( const Array<T> & );
		Array & operator -= ( const Array<T> & );
		Array & operator *= ( const Array<T> & );
		Array & operator /= ( const Array<T> & );

		Array & operator += ( const T & );
		Array & operator -= ( const T & );
		Array & operator *= ( const T & );
		Array & operator /= ( const T & );

		Array   operator + () const;
		Array   operator - () const;

		T & operator [] ( std::size_t );
		T   operator [] ( std::size_t ) const;
		//  TODO rajouter les autres opérateur [] de la classe valarray

		// GETTER
		std::size_t size () const;
		T at ( std::size_t ) const;

		virtual std::size_t getBin ( T )           const;
		virtual T           getX   ( std::size_t ) const;

		// SETTER

		// METHOD
		T min ()                            const;
		//T min ( std::size_t , std::size_t ) const;
		//T min ( T , T )                     const;

		T max ()                            const;
		//T max ( std::size_t , std::size_t ) const;
		//T max ( T , T )                     const;

		T sum ()                            const;
		//T sum ( std::size_t , std::size_t ) const;
		//T sum ( T , T )                     const;

		std::valarray<T> & data ();
		std::valarray<T>   data () const;

		template <typename F>
		Array<T> apply ( F ) const;

		virtual T eval ( T ) const;

	protected:
		virtual void set ( const Array<T> & );

	private:
		std::valarray<T> fdata;
		T fxInit;
		T fxStep;
};

#include "Array.cpp"

#endif
