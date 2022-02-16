CFLAGS ?= -Wall -Wformat-overflow=0 -Wno-unused-variable -Wno-unused-result
CFLAGS += $$(pkg-config --cflags starpumpi-1.3)
LDLIBS += $$(pkg-config --libs starpumpi-1.3)

all: 
	mpic++ $(CFLAGS) simpic.cpp message.cpp $(LDLIBS)

debug: 
	mpic++ -g $(CFLAGS) simpic.cpp message.cpp $(LDLIBS)

clean:
	rm -f a.out *.o *.dat
