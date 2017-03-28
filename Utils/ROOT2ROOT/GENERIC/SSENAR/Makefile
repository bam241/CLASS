CC=g++ -std=c++11 -shared
LIBCINT=-I${CLASS_include}
LIB=`root-config --ldflags` -I${ROOTSYS}/include -I${CLASS_include} -L${CLASS_lib} -lCLASSpkg -fPIC -fopenmp -lgomp -Wunused-result

all: libSsenar.so

Scenar_t_dict.cxx: ZAI.hxx IsotopicVector.hxx Scenar_t.hxx Linkdef.hxx
	rootcint -f $@ -c ${LIBCINT} -p $^

libSsenar.so: Scenar_t_dict.cxx  ZAI.cxx IsotopicVector.cxx Scenar_t.cxx
	${CC} -o $@  ${LIB} $^

clean:
	rm libSsenar.so *dict*
                                                                    
