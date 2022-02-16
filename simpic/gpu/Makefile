CFLAGS ?= -lmpi

all: 
	nvcc simpic.cu message.cpp $(CFLAGS)

debug: 
	nvcc simpic.cu message.cpp $(CFLAGS)

clean:
	rm -f a.out *.o *.dat
