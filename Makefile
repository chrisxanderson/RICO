PKG = rico
VER = 0.1

LIB  = /home/tomf/R/i686-pc-linux-gnu-library/2.15
TGZ  = ./$(PKG)_$(VER).tar.gz
SO   = $(LIB)/$(PKG)/libs/$(PKG).so

RDIR = $(PKG)/R
RR   = $(RDIR)/*.R
SRC  = $(PKG)/src
CPP  = $(SRC)/*.cpp $(SRC)/*.h

all: $(TGZ) $(SO)

$(SO) : $(TGZ)
	R CMD INSTALL $(PKG) -l $(LIB)/

$(TGZ) : $(CPP) $(RR)
	R CMD build $(PKG)


clean:
	@-rm -f $(TGZ) $(SRC)/*.o $(SRC)/*.so $(SRC)/*~

depend:
	@-cd $(SRC) && g++ -MM -I/usr/share/R/include -DNDEBUG -I . *.cpp > DEPS.mk  || exit 1; 
	@-cd $(SRC) && echo "SOURCES = \\"  >  SOURCES.mk
	@-cd $(SRC) && ls *.cpp  | sed -e 's/$$/ \\/' >> SOURCES.mk

.PHONY: all clean depend