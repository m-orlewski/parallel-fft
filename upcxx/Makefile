CC = upcxx
EXEC = fft
RUNNER = upcxx-run
NODEMAKER = /opt/nfs/config/station204_name_list.sh
NODESFILE = nodes
PROCESSES = 8
LIBS = -lm
SRC = $(wildcard *.cpp)
N = 8
DATA_PATH = data/data.txt

all: $(EXEC) run

$(EXEC):
	UPCXX_GASNET_CONDUIT=udp $(CC) -O2 $(SRC) -o $@ -lm

clean:
	rm $(EXEC)

run:
	$(RUNNER) -shared-heap 256M -n $(PROCESSES) $(shell upcxx-nodes nodes) $(EXEC) $(N) $(DATA_PATH)