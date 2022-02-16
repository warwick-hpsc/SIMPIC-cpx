CFLAGS ?= -lmpi -lcusparse
CFLAGS += $$(pkg-config --cflags starpumpi-1.3)
LDLIBS += $$(pkg-config --libs starpumpi-1.3)

all: 
	nvcc $(CFLAGS) simpic.cu message.cpp $(LDLIBS)

debug: 
	nvcc -g $(CFLAGS) simpic.cu message.cpp $(LDLIBS)

clean:
	rm -f a.out *.o *.dat
