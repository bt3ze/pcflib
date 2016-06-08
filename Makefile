all: pcflib.o opdefs.o opflows.o opgen.o threads.o # test pcflib_new.o

pcflib.o: pcflib.cpp pcflib.h opdefs.h opflows.h opgen.h
	gcc -fPIC pcflib.cpp -c -Wall -Werror -g

opdefs.o: opdefs.cpp opdefs.h pcflib.h
	g++ -fPIC opdefs.cpp -c -Wall -Werror -g

opflows.o: opflows.cpp opflows.h opdefs.h pcflib.h
	g++ -fPIC opflows.cpp -c -Wall -Werror -g

opgen.o: opgen.cpp opdefs.h opgen.h pcflib.h
	g++ -fPIC opgen.cpp -c -Wall -Werror -g

threads.o: cthreadpool/thpool.c cthreadpool/thpool.h
	gcc -fPIC cthreadpool/thpool.c -c -Wall -Werror -g -lpthread

#pcflib_new.o: pcflib_new.c pcflib_new.h opdefs_new.c
#	gcc -fPIC pcflib_new.c -c -Wall -Werror -g

#opdefs_new.o: opdefs_new.c opdefs_new.h pcflib_new.h
#	gcc -fPIC opdefs_new.c -c -Wall -Werror -g

#circuitgraph.o: pcflib.h opdefs.h circuitgraph.h circuitgraph.c
#	gcc -fPIC circuitgraph.c -c -Wall -Werror -g

#test: pcflib.o opdefs.o test.c
#	gcc -fPIC -o test test.c pcflib.o opdefs.o -Wall -Werror -g

#cirgen: pcflib.o opdefs.o cirgen.c
#	gcc -fPIC -o cirgen cirgen.c pcflib.o opdefs.o -Wall -Werror -g

clean: 
	rm pcflib.o opdefs.o opflows.o opgen.o cirgen test 
