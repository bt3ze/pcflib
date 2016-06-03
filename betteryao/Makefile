CXX_CFLAGS = -I$(HOME)/lccyao2/pcflib/ -I$(HOME)/include/ -I$(HOME)/local/include -g -L $(HOME)/lib -L$(HOME)/local/lib -Wno-deprecated -D__STDC_LIMIT_MACROS -DFREE_XOR -DRAND_SEED -DGRR  -D_BETTERYAO -DBETTERYAO -O3  -maes -march=native
CXX        = g++
MPI_CXX    = mpicxx
CC = cc
LIBS       =  -lpbc -lgmp -lcrypto -llog4cxx -lpthread #-lgmpxx 
AESNI = -DAESNI #turn this flag on/off for encryption modes
HEADERS    =  Algebra.h Aes.h Bytes.h Env.h NetIO.h Prng.h Commitment.h \
	 #   Circuit.h ClawFree.h  
OBJS       =  Algebra.o Bytes.o Env.o NetIO.o Prng.o Aes.o pcflib.o opdefs.o opflows.o Commitment.o # pcflib_new.o opdefs_new.o
GARBLES    =  garbling.o  GarbledCircuit.o #

#ot extensions
#OTEXTENSION = ot-ot.o otextension.o
#OT-UTILS =  parse_options.o cbitvector.o
#OT-CRYPTO = crypto.o ecc-pk-crypto.o gmp-pk-crypto.o
#OT-ALSZ =  alsz-ot-ext-rec.o alsz-ot-ext-snd.o ot-ext-rec.o ot-ext-snd.o simpleot.o ot-ext.o 
#OTX = otextension.o
#OT-OBJS = cbitvector.o crypto.o ecc-pk-crypto.o gmp-pk-crypto.o

MAINS      = main.cpp YaoBase.cpp Yao.cpp BetterYao4.cpp
#SPLITMAINS = main.cpp YaoBase.cpp Yao.cpp BetterYao4.cpp splityao/BetterYao5.cpp splityao/BetterYaoEvl.cpp splityao/BetterYaoGen.cpp 
MAINS5   = protocol.cpp YaoBase.cpp BetterYao5.cpp
MAINS5HBC = protocol.cpp YaoBase.cpp BetterYao5-hbc.cpp
MACROS    = macros.cpp
#OTE = OTExtension
#OTOT = OTExtension/ot

#OBJECTS_MIRACL= $(OTE)/util/Miracl/*.o
#MIRACL_PATH= -I$(OTE)/util/Miracl
#MIRACL_LIB = -lpthread $(OTE)/util/Miracl/miracl.a -lssl


all : pcflib by5-hbc # otextension.o ot-ot.o ot-crypto.o ot-util.o betteryao5

#betteryao: $(MAINS) $(OBJS) $(OLDGARB)
#	$(MPI_CXX) -o evl -DEVL_CODE $(CXX_CFLAGS) $^ $(LIBS)
#	$(MPI_CXX) -o gen -DGEN_CODE $(CXX_CFLAGS) $^ $(LIBS)

#betteryao5: $(MAINS5) $(OBJS) $(GARBLES) $(MACROS) $(OTX) $(OT-ALSZ) $(OT-UTILS) $(OT-CRYPTO)
#	$(MPI_CXX) -o evl5 -DEVL_CODE -DYAO5 $(CXX_CFLAGS) $^ -I. $(OBJECTS_MIRACL) $(MIRACL_PATH) $(MIRACL_LIB) $(LIBS) 
#	$(MPI_CXX) -o gen5 -DGEN_CODE -DYAO5 $(CXX_CFLAGS) $^ -I. $(OBJECTS_MIRACL) $(MIRACL_PATH) $(MIRACL_LIB) $(LIBS)

#betteryao5: $(MAINS5) $(GARBLES) $(OBJS) $(MACROS)
#	$(MPI_CXX) -o evl5 -DEVL_CODE -DYAO5 $(CXX_CFLAGS) $^ $(LIBS) 
#	$(MPI_CXX) -o gen5 -DGEN_CODE -DYAO5 $(CXX_CFLAGS) $^ $(LIBS)

by5-hbc: $(MAINS5HBC) $(GARBLES) $(OBJS) $(MACROS)
	$(MPI_CXX) -o evl5hbc -DEVL_CODE -DYAO5 $(CXX_CFLAGS) $^ $(LIBS) 
	$(MPI_CXX) -o gen5hbc -DGEN_CODE -DYAO5 $(CXX_CFLAGS) $^ $(LIBS)

pcflib: pcflib.o
	gcc -shared -fPIC pcflib.o -o pcflib.so

pcflib_new: pcflib_new.o
	gcc -shared -fPIC pcflib_new.o -o pcflib_new.so

#debug-test: debug
#	gcc -o ../pcftest ../pcflib.o ../test.c -g -lpthread

#debug: ../pcflib.c ../pcflib.h
#	gcc -Wall -Werror -c -fPIC ../pcflib.c -g

pcflib.o: ../pcflib.cpp ../pcflib.h
	gcc -Wall -Werror -c -fPIC ../pcflib.cpp -g -DBETTERYAO

opdefs.o: ../opdefs.cpp ../pcflib.h ../opdefs.h 
	gcc -Wall -Werror -c -fPIC ../opdefs.cpp -g -DBETTERYAO

opflows.o: ../opflows.cpp ../pcflib.h ../opflows.h ../opdefs.h
	gcc -Wall -Werror -c -fPIC ../opflows.cpp -g -DBETTERYAO

#alsz.o:  $(OTE)/ot/*.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/*.cpp

#ot-util.o: $(OTE)/util/*.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/util/*.cpp

#ot-crypto.o: $(OTE)/util/crypto/*.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/util/crypto/*.cpp

#ot-ot.o: ot-ext.o ot-ext-snd.o ot-ext-rec.o simpleot.o alsz-ot-ext-rec.o alsz-ot-ext-snd.o

#ot-ext.o: $(OTE)/ot/ot-ext.h $(OTE)/ot/ot-ext.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/ot-ext.cpp

#ot-ext-snd.o: $(OTE)/ot/ot-ext.h $(OTE)/ot/ot-ext-snd.h $(OTE)/ot/ot-ext-snd.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/ot-ext-snd.cpp

#ot-ext-rec.o: $(OTE)/ot/ot-ext.h $(OTE)/ot/ot-ext-rec.h $(OTE)/ot/ot-ext-rec.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/ot-ext-rec.cpp

#simpleot.o: $(OTE)/ot/baseOT.h $(OTE)/ot/simpleot.h $(OTE)/ot/simpleot.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/simpleot.cpp

#alsz-ot-ext-rec.o: $(OTE)/ot/ot-ext.h $(OTE)/ot/ot-ext-rec.h $(OTE)/ot/alsz-ot-ext-rec.h $(OTE)/ot/alsz-ot-ext-rec.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/alsz-ot-ext-rec.cpp

#alsz-ot-ext-snd.o: $(OTE)/ot/ot-ext.h $(OTE)/ot/ot-ext-snd.h $(OTE)/ot/alsz-ot-ext-snd.h $(OTE)/ot/alsz-ot-ext-snd.cpp
#	$(CC) -c -IOTExtension/util $(OT-CFLAGS) OTExtension/ot/alsz-ot-ext-snd.cpp

#otextension.o: $(OTOT)/alsz-ot-ext-rec.h $(OTOT)/alsz-ot-ext-snd.h otextension.h otextension.cpp
#	$(CC) -c $(OT-CFLAGS) otextension.cpp

server : ipserver.cpp Bytes.o Env.o NetIO.o
	$(CXX) -pthread -o server $(CXX_CFLAGS) $(CXX_LFLAGS) $^ $(LIBS)

GarbledCircuit.o: Aes.h Bytes.h Env.h Prng.h garbling.h GarbledCircuit.h GarbledCircuit.cpp
	$(MPI_CXX) -msse2 -maes $(AESNI) $(CXX_CFLAGS) -lpthread -c GarbledCircuit.cpp

garbling.o: garbling.h Aes.h Bytes.h garbling.cpp
	$(MPI_CXX) -msse2 -maes $(AESNI) $(CXX_CFLAGS) -lpthread -c garbling.cpp

Env.o : Algebra.h Bytes.h NetIO.h Env.h Env.cpp  #Circuit.h   ClawFree.h 
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
	$(CXX) $(CXX_CFLAGS) -c $(AESNI) Aes.cpp -lm -lrt -lpthread -maes -msse4 -march=native

Bytes.o: Bytes.h Hash.h Bytes.cpp
	$(CXX) $(CXX_CFLAGS) -c Bytes.cpp

Commitment.o: Bytes.h Hash.h Prng.h Commitment.h Commitment.cpp
	$(CXX) $(CXX_CFLAGS) -c Commitment.cpp

clean :
	rm -f *.o gen evl sim test-circuit server
