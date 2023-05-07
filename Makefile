CC = mpicc
EXEC = fft
RUNNER = mpiexec
NODEMAKER = /opt/nfs/config/station204_name_list.sh
NODESFILE = nodes
LIBS = -lm
SRC = $(wildcard *.c)

all: $(EXEC) run

$(EXEC):
	$(CC) $(SRC) -o $@ -lm

clean:
	rm $(EXEC) $(NODESFILE) output/*

run:
	$(NODEMAKER) 1 16 > $(NODESFILE)
	$(RUNNER) -f $(NODESFILE) -n 8 ./$(EXEC) 