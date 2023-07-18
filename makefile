CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type

BIN_DIR = ./bin
BIO_DIR = ./bioUtils
INCLUDES = -I$(BIO_DIR)
BIO_LIBS = -L$(BIO_DIR)/ -lbiotools

all:
	cd $(BIO_DIR); make
	make retroSeeker
	
retroSeeker: retroSeeker.o retroSeekerMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o ${BIN_DIR} retroSeekerMain.o retroSeeker.o \
	$(BIO_LIBS) $(LIBS) 

retroSeeker.o: retroSeeker.cpp retroSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c retroSeeker.cpp
	
retroSeekerMain.o: retroSeekerMain.cpp retroSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c retroSeekerMain.cpp
	
clean:
	rm -f *.o
