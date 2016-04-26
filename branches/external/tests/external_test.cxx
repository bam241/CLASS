#include <gtest/gtest.h>
#include "Array_test.inl"

int main(int argc,char * argv[]) {
 ::testing::InitGoogleTest(&argc,argv);
 return RUN_ALL_TESTS();
}

// COMPILATION
// g++ -o external_test external_test.cxx -lgtest -I $CLASS_include -L $CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp -lTMVA
