all: pcflib.o test #pcflib_new.o

pcflib.o: pcflib.c pcflib.h opdefs.h
	gcc -fPIC pcflib.c -c -Wall -Werror -g

opdefs.o: opdefs.c opdefs.h pcflib.h
	gcc -fPIC opdefs.c -c -Wall -Werror -g

#pcflib_new.o: pcflib_new.c pcflib_new.h opdefs_new.c
#	gcc -fPIC pcflib_new.c -c -Wall -Werror -g

#opdefs_new.o: opdefs_new.c opdefs_new.h pcflib_new.h
#	gcc -fPIC opdefs_new.c -c -Wall -Werror -g

test: pcflib.o opdefs.o test.c
	gcc -fPIC -o test test.c pcflib.o opdefs.o -Wall -Werror -g

cirgen: pcflib.o opdefs.o cirgen.c
	gcc -fPIC -o cirgen cirgen.c pcflib.o opdefs.o -Wall -Werror -g
