SONAME := libmvar.so.0
SHARED_LIB := libmvar.so.0.1.1
STATIC_LIB := libmvar.a

CFLAGS := -Wall -O3 -g -ansi -pedantic -fPIC
OBJS := mvardebug.o mvardie.o mvarsim.o mvarmat.o mvarfit.o mvartest.o

all: $(SHARED_LIB) $(STATIC_LIB) test
.PHONY: all

$(SHARED_LIB): $(OBJS)
	$(CC) -shared -Wl,-soname,$(SONAME) -o $(SHARED_LIB) $(OBJS) -lgsl -lm

$(STATIC_LIB): $(OBJS)
	ar rcs $(STATIC_LIB) $(OBJS)

test: $(STATIC_LIB) test.o
	$(CC) test.o -o test -L. -lmvar -lgsl -lgslcblas
	@LD_LIBRARY_PATH=. ./test

clean:
	rm -rf $(SONAME) $(SHARED_LIB) $(STATIC_LIB)
	rm -rf *.o *.a 
	rm -rf test
.PHONY: clean

deb:
	debuild -i -us -uc -b
.PHONY: deb

install:
	install mvar.h $(DESTDIR)/usr/include
	install $(SHARED_LIB) $(DESTDIR)/usr/lib
	install $(STATIC_LIB) $(DESTDIR)/usr/lib
.PHONY: install

%.o: %.c
	$(CC) $(CFLAGS) -c $<
