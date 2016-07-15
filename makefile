CC=g++
FLAG=-std=c++11
PATH_SRC=./src/
PATH_MATH=./src/math/
PATH_TEST=./test/
RBM=$(PATH_SRC)rbm.h $(PATH_SRC)rbm.cc
MATH=$(PATH_MATH)*

all: rbm_test.o

rbm_test.o: $(MATH) $(RBM) $(PATH_TEST)rbm_test.cc
	$(CC) $(FLAG) $(MATH) $(RBM) $(PATH_TEST)rbm_test.cc -o rbm_test.o