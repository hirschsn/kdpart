
.SUFFIXES:
.SUFFIXES: .cc .o

CXX = mpic++
CXXFLAGS = --std=c++14 -O3 -march=native

AR = ar
RANLIB = ranlib

LIBKDPART_HDR = kdpart.h util/find.h util/codim_sum.h util/mpi_global_vector.h
LIBKDPART_SRC = kdpart.cc util/mpi_global_vector.cc
LIBKDPART_OBJ = $(LIBKDPART_SRC:.cc=.o)
LIBKDPART = libkdpart.a
LIBKDPART_SO = libkdpart.so


TGT_SRC = kdpart_test_par.cc
TGT_OBJ = $(TGT_SRC:.cc=.o)
TGT = kdpart_test_par

all: $(LIBKDPART) $(LIBKDPART_SO) $(TGT)


$(LIBKDPART): $(LIBKDPART_OBJ)
	$(AR) rc $@ $?
	$(RANLIB) $@

$(LIBKDPART_SO): $(LIBKDPART_OBJ)
	$(CXX) -shared -o $@ $?

$(LIBKDPART_OBJ): $(LIBKDPART_HDR)


$(TGT_OBJ): $(LIBKDPART_HDR)

$(TGT): $(LIBKDPART)

.cc.o:
	$(CXX) -fPIC $(CXXFLAGS) -o $@ -c $<

.o:
	$(CXX) -static $(CXXFLAGS) -L. -I. -o $@ $< -lkdpart

clean:
	rm -f $(TGT_OBJ) $(TGT) $(LIBKDPART_SO) $(LIBKDPART) $(LIBKDPART_OBJ)

.PHONY: all clean
