CC=g++

all: sim

sim: main.o cache.o memory.o
	$(CC) -g -o $@ $^

main.o: cache.h

cache.o: cache.h def.h

memory.o: memory.h

.PHONY: clean

clean:
	rm -rf sim *.o
