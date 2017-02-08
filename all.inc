#General make file settings
CC = cc
CPP = c++
LIBTOOL = ar -r

#optimization options
OPTIMOPT=-g
#OPTIMOPT=-O

#C++ 11 support on
STD_OPT = -std=c++14


# Auxilary headers
SNOWGOOSE=$(ROOT)/../snowgoose
COMPI=$(ROOT)/../COMPI
EXTRA_INC = -I$(SNOWGOOSE) -I$(COMPI) 


#common options
COMOPTS = $(STD_OPT)\
	   -I$(ROOT)\
           $(EXTRA_INC)\
           $(ADDIT_INC)\
           $(OPTIMOPT)

#C options
COPTS=$(COMOPTS)  $(PCOPTS)

#C++ options
CPPOPTS= $(COMOPTS) $(ADDIT_INCS) $(PCPPOPTS)

#Linker options for static linking
LDFLAGS = -lm 


MFLAGS =  COPTS="$(COPTS)" CPPOPTS="$(CPPOPTS) LIBTOOL="$(LIBTOOL)"


all: $(BINS) $(LIBS)

tests:
	@for i in $(TESTS); do if ./$$i > /dev/null; then echo TEST PASSED; continue; else echo TEST FAILED; fi done

# Uncomment this if you want testing finishes after first finished test
#	@for i in $(TESTS); do if ./$$i > /dev/null; then echo TEST PASSED; continue; else echo TEST FAILED; return 1; fi done
	

$(BINS): $(LIBS) $(EXTLIBS)

$(LIBS): $(OBJS)
	$(LIBTOOL) $@ $(OBJS)

clean: 
	rm -f *.exe *.o *.a *~ *.log deps.inc

dep:
	$(CC) $(COPTS) -MM -c *.c > deps.inc;\
	$(CPP) $(CPPOPTS) -MM -c *.cpp >> deps.inc;\
        true

.o.exe:
	$(CPP) $(OPTS) -o $@ $< $(LIBS) $(EXTLIBS) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPOPTS) -c $<

.c.o:
	$(CC) $(COPTS) -c $<

.SUFFIXES:
.SUFFIXES: .o .a .cpp .c .exe
