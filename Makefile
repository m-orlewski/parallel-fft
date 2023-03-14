CC=gcc
EXEC = fft
LIBS=-lm
SRC = $(wildcard *.c)

all: $(EXEC) run

$(EXEC):
	$(CC) $(SRC) -o $@ -lm

clean:
	rm fft output/output.txt

run:
	./$(EXEC)