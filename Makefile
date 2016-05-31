all: pcflib.o  # test pcflib_new.o

pcflib.o: pcflib.c pcflib.h opdefs.h
	gcc -fPIC pcflib.c -c -Wall -Werror -g

opdefs.o: opdefs.cpp opdefs.h pcflib.h
	g++ -fPIC opdefs.cpp -c -Wall -Werror -g

#pcflib_new.o: pcflib_new.c pcflib_new.h opdefs_new.c
#	gcc -fPIC pcflib_new.c -c -Wall -Werror -g

#opdefs_new.o: opdefs_new.c opdefs_new.h pcflib_new.h
#	gcc -fPIC opdefs_new.c -c -Wall -Werror -g

circuitgraph.o: pcflib.h opdefs.h circuitgraph.h circuitgraph.c
	gcc -fPIC circuitgraph.c -c -Wall -Werror -g

test: pcflib.o opdefs.o test.c
	gcc -fPIC -o test test.c pcflib.o opdefs.o -Wall -Werror -g

cirgen: pcflib.o opdefs.o cirgen.c
	gcc -fPIC -o cirgen cirgen.c pcflib.o opdefs.o -Wall -Werror -g

clean: 
	rm pcflib.o opdefs.o cirgen test 
