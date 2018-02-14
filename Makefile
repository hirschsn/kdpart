
CC=mpic++
CFLAGS=--std=c++14 -O3 -march=native

INC=kdpart.h util/find.h util/codim_sum.h util/mpi_global_vector.h

TGT=kdpart_test_par

all: $(TGT)

$(TGT): %: %.cc $(INC)
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f $(TGT)

.PHONY: all clean