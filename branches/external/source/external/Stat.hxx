#ifndef CLASS_stat_h
#define CLASS_stat_h

#include <iostream>
#include <cmath>

/*!
	\class Stat<T>
	\brief Class to compute stat on container or on value which arrived progressively

#
La classe permet d'effectuer des statistiques sur un échantillion lors de la
construction ou en ajoutant au fur et à mesure des données.
Cette deuxième possibilité limite les valeurs statistiques calculables :

* Le minimum : $\min_{i} x_i$
* Le maximum : $\max_{i} x_i$
* L'amplitude : $\max_{i}x_i - \min_{i}x_i$
* La moyenne : $\bar{x} = \frac{1}{n}\sum_i x_i$
* L'écart-type (RMS) : $\sigma = \sqrt{ \frac{1}{n} \sum_{i}(x_i - \bar{x})^2 } = \sqrt{ \frac{1}{n}\sum_{i} x_i^2 - \bar{x}^2 }$

Petite précision, RMS signifie *Root Mean Square* ou moyenne quadratique en français. Celle-ci est identique à l'écart-type si la moyenne est nulle puisque la moyenne quadratique se calcule comme suit : $rms = \sqrt{\frac{1}{n}\sum_{i} x_i^2}$
!*/
template <typename T>
class Stat
{
	public:
// CONSTRUCTOR
		Stat ();
		Stat ( const Stat<T> & );
		template <typename ITR>
		Stat ( const ITR , const ITR );

// DESTRUCTOR
		~Stat ();

// OPERATOR
		Stat<T> & operator =  ( const Stat<T> & );
		Stat<T> & operator += ( const T );

// GETTER
		inline std::size_t n () const;

		inline T min () const;
		inline T max () const;
		inline T amp () const;

		inline T avg () const;
		inline T rms () const;

// SETTER
		inline void add ( const T );

	private:
// ATTRIBUTS
		std::size_t fn;    // number of elements
		T           fmin;  // minimum
		T           fmax;  // maximum
		T           fsum;  // sum of all
		T           fsum2; // sum of square
};

#include "stat.cpp"

#endif
