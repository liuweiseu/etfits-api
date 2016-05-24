CC=g++
CFLAGS=-I$(IDIR) -I/usr/include/python2.6 -Wall -fPIC 
LDFLAGS = -shared
LIBS= -lm -lcfitsio -lpython2.6
IDIR=/usr/include/cfitsio

TARGET = test
MODULE = _s6fits.so

$(TARGET):  main.c libs6fits.a  
	$(CC) -L. $(CFLAGS) -I. $^ -o $@ $(LIBS) 

main.o:  main.c  
	$(CC) -L. $(CFLAGS) -I. -c $< -o $@ $(LIBS)  

libs6fits.a: s6fits.o
	ar rc $@ $^
	ranlib libs6fits.a

$(MODULE): s6fits.o s6fits_wrap.o
	$(CC) ${LDFLAGS} -o $@ $^ $(LIBS)

s6fits.o: s6fits.c s6fits.h
	$(CC) ${LDFLAGS} $(CFLAGS) -c -o $@ $< $(LIBS) 

s6fits_wrap.o: s6fits_wrap.cxx 
	$(CC) ${LDFLAGS} $(CFLAGS) -c -o $@ $< $(LIBS) 

clean:
	rm -f *.o *.a *.so test  
