CC=gcc
CFLAGS = -O3 -Wall -Wextra -std=c99 -fopenmp -march=haswell -ffast-math
LDFLAGS=-lpthread -lm

TARGET=pi.out
OBJS=pi.o montecarlo.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

main.o: main.c montecarlo.h
	$(CC) $(CFLAGS) -c pi.c

montecarlo.o: montecarlo.c montecarlo.h
	$(CC) $(CFLAGS) -c montecarlo.c

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
