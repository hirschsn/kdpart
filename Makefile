
.SUFFIXES:
.SUFFIXES: .cc .o

CXX = mpic++
CXXFLAGS = --std=c++14 -O3 -march=native

AR = ar
RANLIB = ranlib

MPIEXEC = mpiexec

LIBKDPART_MASTER_HDR = kdpart.h
LIBKDPART_UTIL_HDR = util/find.h util/codim_sum.h util/mpi_global_vector.h
LIBKDPART_HDR = $(LIBKDPART_MASTER_HDR) $(LIBKDPART_UTIL_HDR)
LIBKDPART_SRC = kdpart.cc util/mpi_global_vector.cc
LIBKDPART_OBJ = $(LIBKDPART_SRC:.cc=.o)
LIBKDPART = libkdpart.a
LIBKDPART_SO = libkdpart.so


TGT_SRC = kdpart_test_par.cc
TGT_OBJ = $(TGT_SRC:.cc=.o)
TGT = kdpart_test_par

all: $(LIBKDPART_SO)

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

check: $(LIBKDPART) $(TGT)
	for i in `seq 1 4`; do $(MPIEXEC) --oversubscribe -n $$i ./kdpart_test_par; done

install: all
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/include/kdpart/util
	cp -f $(LIBKDPART_SO) $(PREFIX)/lib
	cp -f $(LIBKDPART_MASTER_HDR) $(PREFIX)/include/kdpart
	cp -f $(LIBKDPART_UTIL_HDR) $(PREFIX)/include/kdpart/util

.PHONY: all clean check install
