#ifndef _GRAPH_H
#define _GRAPH_H

template <typename T>
class Graph : Array<T>
{
	public:
		// CONSTRUCTOR
		Graph ();
		Graph ( const Graph<T> & );
		// TODO : constructeurs

		// DESTRUCTOR
		virtual ~Graph ();

		// OPERATOR
		virtual Graph<T> & operator = ( const Graph<T> & );

		// GETTER
		virtual std::size_t getBin ( T )           const;
		virtual T           getX   ( std::size_t ) const;

		// METHOD
		T eval ( T ) const;

	protected:
		virtual void set ( const Graph<T> & );

	private:
		std::valarray<T> ftime;
};

#endif
