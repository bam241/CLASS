#ifndef _GRAPH_HXX
#define _GRAPH_HXX

#include <iostream>
#include <algorithm>
#include <iterator>
#include <valarray>

#include <TObject.h>
#include <TGraph.h>


class Graph : public TObject
{
	public:
	// constructor -----------------------------------------------------------
		Graph ();
		Graph ( const Graph & );
		Graph ( std::size_t );
		Graph ( std::size_t , const double * , const double * y=nullptr );
		Graph ( const TGraph & );
		Graph ( const TGraph * );

	// destructor ------------------------------------------------------------
		virtual ~Graph ();

	// operator --------------------------------------------------------------
		Graph &  operator =  ( const Graph & );
		void     operator = ( const TGraph & );
		void     operator = ( const TGraph * );
		operator TGraph () const; // cast operator

		Graph & operator += ( const Graph & );
		Graph & operator -= ( const Graph & );
		Graph & operator *= ( const Graph & );

		const double & operator [] ( std::size_t ) const;
		double & operator [] ( std::size_t );

	// getter ----------------------------------------------------------------
		std::size_t size () const;

		const double * Xbegin () const;
		double *       Xbegin ();
		const double * Ybegin () const;
		double *       Ybegin ();

		const double * Xend   () const;
		double *       Xend   ();
		const double * Yend   () const;
		double *       Yend   ();

		const double * GetX () const;
		double * GetX ();
		const double * GetY () const;
		double * GetY ();

		double GetX ( std::size_t ) const;
		double GetY ( std::size_t ) const;
		double GetY ( double      ) const;

		std::size_t GetN () const;

		int    GetPoint ( int , double & , double & ) const;


	// setter ----------------------------------------------------------------
		void resize ( std::size_t , double x=0 , double y=0 );

		void SetX ( std::size_t , double );
		void SetY ( std::size_t , double );

		void SetPoint ( std::size_t , double , double );

	// method ----------------------------------------------------------------
		double Eval    ( double x ) const;
		double ExpEval ( double x ) const;

		void Sort ();


		ClassDef(Graph,0)
	private:
	// attribut --------------------------------------------------------------
		std::valarray<double> fx;
		std::valarray<double> fy;
};

// external operator ---------------------------------------------------------
Graph operator + ( const Graph & , const Graph & );
Graph operator - ( const Graph & , const Graph & );

Graph operator * ( const Graph & , const Graph & );
Graph operator * ( double , const Graph & );
Graph operator * ( const Graph & , double );

Graph operator / ( const Graph & , double );

#endif
