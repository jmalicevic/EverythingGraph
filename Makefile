#CC=clang-3.5
CFLAGS=-fcilkplus -lcilkrts -DCILK -pthread  -O2 -std=gnu99 -fms-extensions -ggdb3 
#CFLAGS=-O2 -fopenmp -pthread -std=gnu99 -fms-extensions -ggdb3
CC=g++
LDFLAGS= -lm  -llapack -lblas -lnuma

OBJ= random.o utils.o buffer.o  init_all.o
ALGOS=sssp_pushpull spmv bfsgrid_cilk wcc prgrid_cilk pagerank_simple bfs_simple bfs_numa pr_numa 
DIRS= $(patsubst %, .%, $(ALGOS))
MK= $(patsubst %, .%/Makefile, $(ALGOS))

DIRS= $(patsubst %, .%, $(ALGOS))
MK= $(patsubst %, .%/Makefile, $(ALGOS))
export

.PHONY:clean all $(ALGOS)

all: $(ALGOS) 

$(MK):Makefile.tmpl
	@(echo "#Updating makefile $@"; mkdir -p $(dir $@); cp Makefile.tmpl $@;)

$(ALGOS): %: .%/Makefile
	@(make -C .$@)

clean:
	rm -rf *.o $(DIRS) $(ALGOS)

