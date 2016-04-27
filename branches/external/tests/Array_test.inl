#include <iostream>
#include <algorithm>
#include "external/Array.hxx"
#include "external/Graph.hxx"

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

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

TEST ( TestArray , operatorEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a1(n,f);
	const Array<double> a2(2*n,f); // a2 est plus grand et ne possède pas les mêmes valeurs

	a1 = a2;
	EXPECT_EQ ( a1.size() , a2.size() );
	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a1[i] , a2[i] );
	}
}

TEST ( TestArray , operatorPEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a1(n,f);
	const Array<double> a2(n,f); // a2 ne possède pas les mêmes valeurs

	const Array<double> tmp = a1;

	a1 += a2;

	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a1[i] , tmp[i]+a2[i] );
	}
}

TEST ( TestArray , operatorMEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a1(n,f);
	const Array<double> a2(n,f); // a2 ne possède pas les mêmes valeurs

	const Array<double> tmp = a1;

	a1 -= a2;

	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a1[i] , tmp[i]-a2[i] );
	}
}

TEST ( TestArray , operatorSEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a1(n,f);
	const Array<double> a2(n,f); // a2 ne possède pas les mêmes valeurs

	const Array<double> tmp = a1;

	a1 *= a2;

	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a1[i] , tmp[i]*a2[i] );
	}
}

TEST ( TestArray , operatorDEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=1.5; f.h = 1;

	Array<double> a1(n,f);
	const Array<double> a2(n,f); // a2 ne possède pas les mêmes valeurs

	const Array<double> tmp = a1;

	a1 /= a2;

	EXPECT_EQ ( a1.size() , a2.size() );
	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a1[i] , tmp[i]/a2[i] );
	}
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , operatorVPEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a(n,f);
	const double v = 3.141592653589;

	const Array<double> tmp = a;

	a += v;

	for ( unsigned int i = 0 ; i < a.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a[i] , tmp[i]+v );
	}
}

TEST ( TestArray , operatorVMEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a(n,f);
	const double v = 3.141592653589;

	const Array<double> tmp = a;

	a -= v;

	for ( unsigned int i = 0 ; i < a.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a[i] , tmp[i]-v );
	}
}

TEST ( TestArray , operatorVSEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a(n,f);
	const double v = 3.141592653589;

	const Array<double> tmp = a;

	a *= v;

	for ( unsigned int i = 0 ; i < a.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a[i] , tmp[i]*v );
	}
}

TEST ( TestArray , operatorVDEQ )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=1; f.h = 1;

	Array<double> a1(n,f);
	double v = 3.141592653589;

	Array<double> tmp = a1;

	a1 /= v;

	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a1[i] , tmp[i]/v );
	}
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , operatorPLUS )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a1(n,f);
	Array<double> a2;

	a2 = +a1;

	EXPECT_EQ ( a2.size() , a1.size() );
	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a2[i] , +a1[i] );
	}
}

TEST ( TestArray , operatorMINUS )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a1(n,f);
	Array<double> a2;

	a2 = -a1;

	EXPECT_EQ ( a2.size() , a1.size() );
	for ( unsigned int i = 0 ; i < a1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( a2[i] , -a1[i] );
	}
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , operatorCROCHET )
{
	// on a testé les crochets en lecture mais pas en écriture
	const std::size_t n = 1<<5;
	const double e = 2.71828182;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	Array<double> a(n,f);

	a[5] = e;

	EXPECT_DOUBLE_EQ ( a[5] , e );
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , size )
{
	const std::size_t n = 1<<5;

	const Array<double> a1;
	const Array<double> a2(n);

	EXPECT_EQ ( a1.size() , 0 );
	EXPECT_EQ ( a2.size() , n );
}
TEST ( TestArray , at )
{
	const std::size_t n = 1<<5;

	Array<unsigned int> a(n);

	for ( unsigned int i = 0 ; i<n ; ++i )
		{ a[i] = i; }

	EXPECT_EQ ( a.at(5) , 5 );
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , getBinX )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	const Array<double> a( n , f , -1 , 0.01 );
	//               ^   ^    ^    ^
	//            size fonct xInit xStep

	EXPECT_EQ ( a.getBin(-0.9) , 9 );
	EXPECT_EQ ( a.getBin(-1)   , 0 );
	EXPECT_DOUBLE_EQ ( a.getX(10) , -0.9 );
	EXPECT_DOUBLE_EQ ( a.getX(0) , -1   );
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , minmax )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=(int)-n/2; f.h = 1;

	const Array<double> a( n , f , -1 , 0.01 );

	EXPECT_DOUBLE_EQ ( a.min() , a[0]   );
	EXPECT_DOUBLE_EQ ( a.max() , a[n-1] );
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , sum )
{
	const std::size_t n = 1<<5;
	fonctor f;
	f.value=0; f.h = 1;

	const Array<double> a( n , f , -1 , 0.01 );

	EXPECT_DOUBLE_EQ ( a.sum() , n*(n+1)/2 );
}
//////////////////////////////////////////////////////////////////////////////
TEST ( TestArray , operators )
{ // ceci n'est pas un test rigoureux
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

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

TEST ( TestGraph , constructor )
{
	const std::size_t n = 1<<5;

	std::vector<double> v(n);
	fonctor f;
	std::generate(v.begin(),v.end(),f);

	Graph<double> g1(n,&v[0]);
	Graph<double> g2 = g1;

	Graph<double> g3 = g1 + g2;

	EXPECT_EQ ( g1.size() , n );
	
	for ( unsigned int i = 0 ; i < n ; ++i )
	{
		EXPECT_DOUBLE_EQ ( g3[i] , g1[i]+g2[i] );
	}

	g2 += 0.01;
	Graph<bool> tmp = g1 < g2;
	for ( bool * it = tmp.begin() ; it  != tmp.end() ; ++it )
		{ EXPECT_TRUE ( *it ); }
}

//////////////////////////////////////////////////////////////////////////////
TEST ( TestGraph , operatorEQ )
{
	const std::size_t n = 1<<5;

	std::vector<double> time(n);
	fonctor f;
	std::generate(time.begin(),time.end(),f);

	Graph<double> g1(n,f,time);
	const Graph<double> g2(2*n,f,time); // g2 est plus grand et ne possède pas les mêmes valeurs

	g1 = g2;
	EXPECT_EQ ( g1.size() , g2.size() );
	for ( unsigned int i = 0 ; i < g1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( g1[i] , g2[i] );
		EXPECT_DOUBLE_EQ ( g1.getX(i) , g2.getX(i) );
	}
}

TEST ( TestGraph , operatorPEQ )
{
	const std::size_t n = 1<<5;

	std::vector<double> time(n);
	fonctor f;
	f.value=(int)-n/2; f.h = 1;
	std::generate(time.begin(),time.end(),f);

	Graph<double> g1(n,f,time);
	const Graph<double> g2(n,f,time);

	const Array<double> tmp = g1;

	g1 += g2;

	for ( unsigned int i = 0 ; i < g1.size() ; ++i )
	{
		EXPECT_DOUBLE_EQ ( g1[i] , tmp[i]+g2[i] );
	}
}

TEST ( TestGraph , getBinX )
{
	const std::size_t n = 1<<5;

	std::vector<double> v(n);
	fonctor f;
	f.value=2; f.h = 0.1; // x need to be positive
	std::generate(v.begin(),v.end(),f);

	const Graph<double> g(n,f,v);

	const unsigned int i = 5;
	const double x = v[i] + 0.25*(v[i]-v[i+1])/2;

	EXPECT_EQ ( g.getBin(x) , i );
	EXPECT_DOUBLE_EQ ( g.getX(i) , v[i] );
}


