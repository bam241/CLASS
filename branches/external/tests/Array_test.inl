#include "external/Array.hxx"

struct fonctor
{
	int value;
	double h;
	fonctor () :
		value(-1) , h(0.01)
	{ ; }
	
	fonctor ( const fonctor & a ) :
		value(a.value) , h(a.h)
	{ ; }

	~fonctor ()
	{ ; }

	double operator () ()
	{ return (++value)*h; }
};

TEST ( TestArray , constructor )
{
	const std::size_t n = 1<<5;

	Array<int> a1;
	Array<double> a2 ( n , 3.14 , 0.01 );
	Array<double> a3 ( a2 );
	
	// génération d'un tableau et vecteur de référence
	double * tab = new double[n];
	std::vector<double> v(n);
	for ( std::size_t i = 0 ; i < n ; ++i )
	{
		v[i] = tab[i] = 3.141592653589*i;
	}

	Array<double> a4 ( tab , n , -1 , 0.01 );   // constructeur à partir d'un tableau
	Array<double> a5 ( v , -1 , 0.01 );         // constructeur à partir d'un vector<double>
	Array<double> a6 ( a5.data() , -1 , 0.01 ); // constructeur à partir d'un valarray<double>

	fonctor f;
	Array<double> a7 ( n , f , 0 , 0.01 );      // constructeur à partir d'une fonction
	

	EXPECT_EQ ( a1.size() , 0 );
	EXPECT_EQ ( a2.size() , n );
	EXPECT_EQ ( a3.size() , n );
	EXPECT_EQ ( a4.size() , n );
	EXPECT_EQ ( a5.size() , n );
	EXPECT_EQ ( a6.size() , n );
	EXPECT_EQ ( a7.size() , n );
	EXPECT_DOUBLE_EQ ( a7.sum() , 4.96 ); // l'égalité stricte entre deux doubles est toujours foireuse

	// vérification de l'égalité a4 et a5
	Array<bool> tmp = (a4 == a5);
	EXPECT_DOUBLE_EQ ( a4.xBegin() , a5.xBegin() );
	EXPECT_TRUE ( tmp.sum() );

	delete[] tab; tab=0x0;
}

TEST ( TestArray , operator )
{
	const std::size_t n = 1<<5;

	fonctor f;
	f.value=(int)-n/2; f.h = 1;
	Array<int> a1(n,f);
	Array<int> a2 = -a1;

	Array<bool> tmp = (a1*a1) == (a2*a2);

	EXPECT_EQ ( (a2+a1).sum() , 0 );
	EXPECT_TRUE ( tmp.sum() );
	EXPECT_LT ( *(a1.begin()) , *(a1.end()) );

	for ( unsigned int i = 0 ; i < n ; ++i )
	{ // test all value
		EXPECT_EQ ( a1[i] , -a2[i] );
	}

	tmp = (a1+a1) == (2*a1);
	for ( unsigned int i = 0 ; i < n ; ++i )
	{
		EXPECT_TRUE ( tmp[i] );
	}

	tmp = (a1+a1) == (a1*2);
	for ( unsigned int i = 0 ; i < n ; ++i )
	{
		EXPECT_TRUE ( tmp[i] );
	}
}

TEST ( TestArray , minmaxsum )
{
	const std::size_t n = 1<<5;

	fonctor f;
	f.value=(int)-n/2; f.h = 1;
	Array<int> a1(n,f);
	Array<int> a2 = -a1;

	EXPECT_EQ ( a1.min() , -a2.max() );
	EXPECT_EQ ( a1[0] , a1.min() );

	int s = 0;
	for ( int * i = a1.begin() ; i != a1.end() ; ++i )
		{ s += *i; } // calcul de la somme

	EXPECT_EQ ( a1.sum() , s );
}
