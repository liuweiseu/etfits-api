CC=g++
CFLAGS=-I$(IDIR) 
LIBS= -lm -lcfitsio
IDIR=/usr/include/cfitsio

TARGET = test

$(TARGET):  main.c libs6fits.a 
	$(CC) -L. $(CFLAGS) -I. $^ -o $@ $(LIBS)  

main.o:  main.c  
	$(CC) -L. $(CFLAGS) -I. -c $< -o $@ $(LIBS)  

libs6fits.a: s6fits.o
	ar rc $@ $^
	ranlib libs6fits.a

s6fits.o: s6fits.c s6fits.h
	$(CC) $(CFLAGS) -c -o $@ $< $(LIBS) 

clean:
	rm -f *.o *.a test  
