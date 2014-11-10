# ===========================================================================
#  Makefile free-streming                              Jia Liu Apr. 04, 2014
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := $(shell ../crank/configure_compiler_C++.sh)
CFLAGS= $(shell ../crank/configure_compiler_C++_Flags.sh)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	lm.e
endif

SRC		=	main.cpp arsenal.cpp CellData.cpp EOS.cpp Freestreaming.cpp \
			gauss_quadrature.cpp LdMatching.cpp mistools.cpp ParameterReader.cpp \

INC		= 	arsenal.h CellData.h EOS.h Freestreaming.h \
			gauss_quadrature.h LdMatching.h mistools.h ParameterReader.h\


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(LDFLAGS) $(OBJECTS) -o $(TARGET)
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
./EOS.cpp: arsenal.h
./Table.cpp: arsenal.h
./Freestreaming.cpp: gauss_quadrature.h mistools.h arsenal.h
./LdMatching.cpp: gauss_quadrature.h mistools.h arsenal.h	
./main.cpp: LdMatching.h Freestreaming.h gauss_quadrature.h ParameterReader.h

