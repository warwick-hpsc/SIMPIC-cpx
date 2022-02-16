CFLAGS ?= -lmpi -lcusparse

all: 
	nvcc simpic.cu message.cpp $(CFLAGS)

debug: 
	nvcc -g -G simpic.cu message.cpp $(CFLAGS)

clean:
	rm -f a.out *.o *.dat *.png
