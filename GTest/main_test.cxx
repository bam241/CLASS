#include <gtest/gtest.h>
#include "test_iv.inl"

int main(int argc,char * argv[]) {
	::testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}

// COMPILATION
// g++ -o CLASS_test main_test.cxx -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -lTMVA -lgtest
