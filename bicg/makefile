.PHONY: all debug
all: a.out
CFLAGS= -isystem /usr/lib/x86_64-linux-gnu/openmpi/include -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wcast-align  -Wfloat-equal \
-Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual \
-Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format



%.o: %.cpp
	mpicxx -c -g $(CFLAGS) $< -o $@

main.o: main.cpp 

a.out: main.o
	mpicxx -o a.out main.o 

clean:
	rm -f *.out *.o
#-Werror -Wunused -pedantic -pedantic-errors 
