#include <gtest/gtest.h>
#include "IsotopicVector/test_iv.inl"
#include "Reactor/test_Reactor1.inl"
#include "Fleet/test_PWR_UOX_MOX.inl"

int main(int argc,char * argv[]) {
	testing::InitGoogleTest(&argc,argv);
	return RUN_ALL_TESTS();
}

/*

#### CLEAN

rm RunTest gtest-all.o lib/libgtest.a DecayDataBank.log

#### Google test library compilation

g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR}/ -pthread -c ${GTEST_DIR}/src/gtest-all.cc
ar -rv lib/libgtest.a gtest-all.o

#### Main Test compilation

g++ -isystem ${GTEST_DIR}/include -pthread main_test.cxx lib/libgtest.a -o RunTest -I${GTEST_DIR}/include -I$CLASS_include -I$CLASS_external -I$CLASS_Equivalence -I$CLASS_Irradiation -I$CLASS_XS -L$CLASS_lib -lCLASSpkg `root-config --cflags` `root-config --libs`

#### Run the test 

./RunTest

*/
