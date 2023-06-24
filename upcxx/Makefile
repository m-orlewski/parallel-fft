CC = upcxx
EXEC = fft
RUNNER = upcxx-run
NODEMAKER = /opt/nfs/config/station204_name_list.sh
NODESFILE = nodes
PROCESSES = 9
LIBS = -lm
SRC = $(wildcard *.cpp)

all: $(EXEC) run

$(EXEC):
	UPCXX_GASNET_CONDUIT=udp $(CC) -O2 $(SRC) -o $@ -lm

clean:
	rm $(EXEC) $(NODESFILE)

run:
	$(NODEMAKER) 1 16 > $(NODESFILE)
	$(RUNNER) -shared-heap 256M -n $(PROCESSES) $(shell upcxx-nodes nodes) $(EXEC) 