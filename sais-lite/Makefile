# Makefile for suftest and test

# options
CC						= gcc
#CXX						= g++
#OUTPUT_OPTION	= -o $@
CFLAGS				= -O3 -fomit-frame-pointer -fPIC
#CXXFLAGS			= -O3 -fomit-frame-pointer
CPPFLAGS			= -Wall -DNDEBUG
LDFLAGS				= 
LDLIBS				= 
#TARGET_ARCH		=

LIBPATH = ../lib

# targets
.PHONY: all
all: directory suftest dll
suftest: sais.o suftest.o
dll: suftest
	$(CC) -shared -Wl,-soname,libsais.so -o $(LIBPATH)/libsais.so *.o
test:
	$(CC) -O -g -Wall test.c sais.c -o test
	./test
	$(RM) test test.exe
directory:
	test -d $(LIBPATH) || mkdir $(LIBPATH)
distclean: clean
clean:
	$(RM) suftest suftest.exe test test.exe sais.o suftest.o

# dependencies
sais.o suftest.o: sais.h Makefile
