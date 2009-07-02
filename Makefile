LIBRARY := libmvar.a
CFLAGS := -Wall -O3 -g -ansi -pedantic
OBJS := mvardebug.o mvardie.o mvarsim.o

all: $(LIBRARY) test
.PHONY: all

$(LIBRARY): $(OBJS)
	ar rcs $(LIBRARY) $(OBJS)

test: $(LIBRARY) test.o
	$(CC) test.o -o test -L. -lmvar -lgsl -lgslcblas
	./test

clean:
	rm -rf *.o *.a test
.PHONY: clean

%.o: %.c
	$(CC) $(CFLAGS) -c $<
