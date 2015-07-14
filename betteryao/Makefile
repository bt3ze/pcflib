CXX_CFLAGS = -I$(HOME)/lccyao2/pcflib/ -I$(HOME)/include/ -I$(HOME)/local/include -g -L $(HOME)/lib -L$(HOME)/local/lib -Wno-deprecated -D__STDC_LIMIT_MACROS -DFREE_XOR -DRAND_SEED -DGRR  -D_BETTERYAO -DBETTERYAO -O0
CXX        = g++
MPI_CXX    = mpicxx

LIBS       = -lpbc -lgmp -lcrypto -llog4cxx -lpthread -lgmpxx
HEADERS    = Algebra.h Aes.h Bytes.h Env.h NetIO.h Prng.h Commitment.h \
	 # Circuit.h 	GarbledBase.h GarbledHBC.h GarbledMal.h   ClawFree.h  garbled_circuit.h garbled_circuit_hbc.h garbled_circuit_m.h
OBJS       = Algebra.o Bytes.o Env.o NetIO.o Prng.o Aes.o pcflib.o opdefs.o Commitment.o
OLDGARB	   = garbled_circuit.o garbled_circuit_hbc.o garbled_circuit_m.o  # Circuit.o GarbledCct3.o
GARBLES    = GarbledCircuit.o #GarbledBase.o GarbledHBC.o GarbledMal.o 
OTEXTENSION = alsz.o ot-crypto.o ot-util.o
OT-UTILS =  parse_options.o cbitvector.o
OT-CRYPTO = crypto.o ecc-pk-crypto.o gmp-pk-crypto.o
OT-ALSZ = alsz-ot-ext-rec.o alsz-ot-ext-snd.o ot-ext.o ot-ext-snd.o ot-ext-rec.o simpleot.o
#OT-OBJS = cbitvector.o crypto.o ecc-pk-crypto.o gmp-pk-crypto.o
MAINS      = main.cpp YaoBase.cpp Yao.cpp BetterYao4.cpp
SPLITMAINS = main.cpp YaoBase.cpp Yao.cpp BetterYao4.cpp splityao/BetterYao5.cpp splityao/BetterYaoEvl.cpp splityao/BetterYaoGen.cpp 
MAINS5   = protocol.cpp YaoBase.cpp BetterYao5.cpp
MACROS    = macros.cpp
OTE = OTExtension
#OT-UTIL = OTExtension/util
#OT-CRYPTO = OTExtension/util/crypto
#OT-OT = OTExtension/ot
OBJECTS_MIRACL= $(OTE)/util/Miracl/*.o
MIRACL_PATH= -I$(OTE)/util/Miracl
MIRACL_LIB = -lpthread $(OTE)/util/Miracl/miracl.a -lssl -lcrypto -lgmp -lgmpxx


all : pcflib alsz.o ot-crypto.o ot-util.o betteryao5

betteryao: $(MAINS) $(OBJS) $(OLDGARB)
	$(MPI_CXX) -o evl -DEVL_CODE $(CXX_CFLAGS) $^ $(LIBS)
	$(MPI_CXX) -o gen -DGEN_CODE $(CXX_CFLAGS) $^ $(LIBS)

betteryao5: $(MAINS5) $(OBJS) $(GARBLES) $(MACROS)
	$(MPI_CXX) -o evl5 -DEVL_CODE -DYAO5 $(CXX_CFLAGS) $^ $(OT-UTILS) $(OT-CRYPTO) $(OT-ALSZ) $(OBJECTS_MIRACL) $(MIRACL_PATH) $(MIRACL_LIB) $(LIBS)
	$(MPI_CXX) -o gen5 -DGEN_CODE -DYAO5 $(CXX_CFLAGS) $^ $(OT-UTILS) $(OT-CRYPTO) $(OT-ALSZ) $(OBJECTS_MIRACL) $(MIRACL_PATH) $(MIRACL_LIB) $(LIBS)

pcflib: pcflib.o
	gcc -shared -fPIC pcflib.o -o pcflib.so

#debug-test: debug
#	gcc -o ../pcftest ../pcflib.o ../test.c -g -lpthread

#debug: ../pcflib.c ../pcflib.h
#	gcc -Wall -Werror -c -fPIC ../pcflib.c -g

pcflib.o: ../pcflib.c ../pcflib.h
	gcc -Wall -Werror -c -fPIC ../pcflib.c -g -DBETTERYAO

opdefs.o: ../opdefs.c ../pcflib.h ../opdefs.h 
	gcc -Wall -Werror -c -fPIC ../opdefs.c -g -DBETTERYAO

#sim: $(MAINS) $(OBJS)
#	$(MPI_CXX) -o sim $(CXX_CFLAGS) $^ $(LIBS)

#test-circuit: test-circuit.cpp $(OBJS)
#	$(MPI_CXX) $(CXX_CFLAGS) -o test-circuit $(CXX_LFLAGS) $^ $(LIBS)

alsz.o:  $(OTE)/util/*.cpp
	$(CC) -c -IOTExtension/ot $(OT-CFLAGS) OTExtension/ot/*.cpp

ot-util.o: $(OTE)/util/*.cpp
	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/util/*.cpp

ot-crypto.o: $(OTE)/util/crypto/*.cpp
	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/util/crypto/*.cpp

#parse_options.o: $(OT-UTIL)/parse_options.h  $(OT-UTIL)/parse_options.cpp
#	$(CC) -c -I$(OT-UTIL) $(OT_CFLAGS) $(OT-UTIL)/parse_options.cpp

#cbitvector.o: $(OT-UTIL)/cbitvector.h  $(OT-UTIL)/cbitvector.cpp $(OT-UTIL)/typedefs.h $(OT-CRYPTO)/crypto.h
#	$(CC) -c -I$(OT-UTIL) $(OT_CFLAGS) $(OT-UTIL)/cbitvector.cpp

server : ipserver.cpp Bytes.o Env.o NetIO.o
	$(CXX) -pthread -o server $(CXX_CFLAGS) $(CXX_LFLAGS) $^ $(LIBS)

GarbledCircuit.o: Algebra.o Aes.h Bytes.h Env.h Prng.h GarbledCircuit.h GarbledCircuit.cpp
	$(MPI_CXX) -msse2 $(CXX_CFLAGS) -c GarbledCircuit.cpp

Env.o : Algebra.h Bytes.h ClawFree.h NetIO.h Env.h Env.cpp #Circuit.h 
	$(CXX) $(CXX_CFLAGS) -c Env.cpp

NetIO.o : Bytes.h NetIO.h NetIO.cpp
	$(CXX) $(CXX_CFLAGS) -c NetIO.cpp

Algebra.o: Bytes.h Prng.h Algebra.h Algebra.cpp
	$(CXX) $(CXX_CFLAGS) -c Algebra.cpp

#Circuit.o : Bytes.h Circuit.h Circuit.cpp
#	$(CXX) $(CXX_CFLAGS) -c Circuit.cpp

Prng.o: Bytes.h Prng.h Prng.cpp
	$(CXX) $(CXX_CFLAGS) -c Prng.cpp

Aes.o: Bytes.h Aes.h Aes.cpp
	$(CXX) $(CXX_CFLAGS) -c Aes.cpp

Bytes.o: Bytes.h Hash.h Bytes.cpp
	$(CXX) $(CXX_CFLAGS) -c Bytes.cpp

Commitment.o: Bytes.h Hash.h Prng.h Commitment.h Commitment.cpp
	$(CXX) $(CXX_CFLAGS) -c Commitment.cpp

clean :
	rm -f *.o gen evl sim test-circuit server
