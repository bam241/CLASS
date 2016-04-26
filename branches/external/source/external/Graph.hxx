#ifndef _GRAPH_H
#define _GRAPH_H

template <typename T>
class Graph : Array<T>
{
	public:
		// CONSTRUCTOR
		Graph ();
		Graph ( const Graph<T> & );

		Graph ( std::size_t );
		
		Graph ( std::size_t , const T * );
		Graph ( const T * , std::size_t  , const T * );
		Graph ( const std::valarray<T> & , const T * );
		Graph ( const std::vector<T>   & , const T * );
		template <typename F>
		Graph ( std::size_t , F , const T * );

		// DESTRUCTOR
		virtual ~Graph ();

		// OPERATOR
		virtual Graph<T> & operator = ( const Graph<T> & );

		Graph<T> & operator += ( const Graph<T> & );
		Graph<T> & operator -= ( const Graph<T> & );
		Graph<T> & operator *= ( const Graph<T> & );
		Graph<T> & operator /= ( const Graph<T> & );

		Graph<T> & operator += ( const T & );
		Graph<T> & operator -= ( const T & );
		Graph<T> & operator *= ( const T & );
		Graph<T> & operator /= ( const T & );

		Graph<T>   operator + () const;
		Graph<T>   operator - () const;

		// GETTER
		virtual std::size_t getBin ( T )           const;
		virtual T           getX   ( std::size_t ) const;

		virtual T xBegin () const;
		virtual T xEnd   () const;

		// METHOD
		std::valarray<T> & time ();
		std::valarray<T>   time () const;

		template <typename F>
		virtual Graph<T> apply ( F ) const;

		T eval ( T ) const;

	protected:
		virtual void set ( const Graph<T> & );

	private:
		std::valarray<T> ftime;

};

#endif
