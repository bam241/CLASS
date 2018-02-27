#include <gtest/gtest.h>
#include "test_iv.inl"

int main(int argc,char * argv[]) {
	::testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}

/*
Google test library compilation

g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv libgtest.a gtest-all.o

Main Test compilation
g++ -isystem ${GTEST_DIR}/include -pthread main_test.cxx libgtest.a -o MyTest -I${GTEST_DIR}/include -I$CLASS_include -I$CLASS_external -I$CLASS_Equivalence -I$CLASS_Irradiation -I$CLASS_XS -L$CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs`
*/
