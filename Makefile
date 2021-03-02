##################################################################
#
#   makefile
# 
#   for the MARS software
#
##################################################################
# @maintitle

# @code

#
#  please change all system depend values in the 
#  config.mk.${OSTYPE} file 
#
#
include Makefile.conf.$(OSTYPE)
include Makefile.conf.general

#
# PROGRAMS = readraw merpp mars test mona status
# mona now disabled (see Changelog) 
#
PROGRAMS1 = ape callisto caspar celestina coach datacubes electronflux fakerecalib fitebl flute flute_mod foam fold fluxlc inforc iscream  matelsim made-up mars \
	   melibea merpp mola odie online osteria pasta printhardware psearch quate readdaq \
	   readraw selectmc showlog showplot sinope sorcerer star superstar truee zinc 

ifdef MDMSYS
    PROGRAMS  = $(PROGRAMS1) glikeInputs
else
    PROGRAMS = $(PROGRAMS1)
endif

ifneq ($(OSTYPE),darwin)
    SOLIB    = libmars.so
endif

CINT     = M

#
#  connect the include files defined in the config.mk file
#
#    WARNING: the result (whether the linkage works or not) depends on the
#             order of the libraries. It seems, that the most base library
#             must be the last one
#

#
#  ----->>>   mars libraries
#
SUBDIRS1 = manalysis \
          mastro \
          matca \
          mbadpixels \
	  mbase \
	  mcalib \
          mcamera \
	  mcta \
          mdata \
          mdatacheck \
          mdisp \
          mfbase \
          mfileio \
          mfilter \
          mflux \
          mgeom \
          mgui \
          mhbase \
          mhcalib \
          mhflux \
          mhft \
          mhist \
          mhistmc \
          mhvstime \
          mimage \
	  mjobs \
	  mmain \
          mmc \
	  mmctelsim \
          mmontecarlo \
          mmuon \
	  monline \
          mpedestal \
          mpointing \
          mpulsar   \
          mranforest \
          mraw \
          mreport \
          msignal \
          mskyplot \
          msorcerer \
          msql \
	  mstarcam \
          msupercuts \
          mtools \
	  mtrigger \
	  mtruee

ifdef MDMSYS
#    MDMLIB    = $(MDMSYS)/libmdm.so
    SUBDIRS=$(SUBDIRS1) $(MDMSYS)/source
else
    SUBDIRS=$(SUBDIRS1)
endif


#LIBRARIES = $(SUBDIRS:%=lib/lib%.a)
LIBRARIES = $(SUBDIRS:=.a)
MRPROPERS = $(SUBDIRS:=.mrproper)
CLEANERS  = $(SUBDIRS:=.clean)
ifneq ($(OSTYPE),darwin)
LIBS      = $(SOLIB)
else
LIBS      = $(DYLIB)
endif

#------------------------------------------------------------------------------
.SUFFIXES: .c .cc .h .o 

SRCFILES = 

############################################################
all: version $(LIBS) $(PROGRAMS)
	@echo " Done. "
	@echo " "

version: 
	@rm -f mbase/MMarsVersion.o

static: LIBS=$(SUBDIRS:=/*.o) $(OBJS)
#static: rmlib $(LIBRARIES) $(PROGRAMS)
static: $(LIBRARIES) $(PROGRAMS)
	@echo " Done. "
	@echo " "

include Makefile.rules

#
# Use $(CXX) -v ... for a more verbose output
#

# This is a special workaround to create the shared object (bundle, plugin)
# for root and the dynlib (to be linked with the executable) on Mac OSX
ifneq ($(OSTYPE),darwin)

# ROOTGLIBS must be there - why? How can I link the libraries?
$(SOLIB): $(LIBRARIES) $(OBJS) $(HEADERS)
	@echo " Linking shared object $(SOLIB) ..."
	$(CXX) $(CXXFLAGS) $(SOFLAG) $(OBJS) $(SUBDIRS1:=/*.o) $(ROOTGLIBS) -o $@

$(PROGRAMS): $(SOLIB) $(PROGRAMS:=.o)
	@echo " Linking $@ ..." 
	$(CXX) $(CXXFLAGS) $@.o $(MARS_LIB) -o $@ $(ROOTGLIBS) $(SOLIB) $(MDMLIB)

#glikeInputs: $(SOLIB) $(MDMLIB) glikeInputs.o
#	@echo " Linking glikeInputs ..." 
#	$(CXX) $(CXXFLAGS) glikeInputs.o $(MARS_LIB) -o $@ $(ROOTGLIBS) $(SOLIB) $(MDMLIB)

# Use this to link the programs statically - for gprof
#$(PROGRAMS): $(OBJS) $(HEADERS) $(PROGRAMS:=.o)
#	@echo " Linking $@ ..." 
#	$(CXX) $(CXXFLAGS) $(ROOTGLIBS) $(OBJS) $(SUBDIRS:=/*.o) $@.o $(MARS_LIB) -o $@
else
$(DYLIB): $(LIBRARIES) $(OBJS) $(HEADERS)
	@echo " Linking dylib $(DYLIB) ..."
	$(CXX) $(CXXFLAGS) $(DYFLAG) $(OBJS) $(SUBDIRS1:=/*.o) $(ROOTGLIBS) -o $@

$(PROGRAMS): $(DYLIB) $(PROGRAMS:=.o)
	@echo " Linking mac executable $@ ..." 
	$(CXX) $(CXXFLAGS) -bind_at_load $(ROOTGLIBS) $(DYLIB) $@.o $(MARS_LIB) $(MDMLIB) -o $@
endif


ifneq ($(OSTYPE),darwin)
 dox: $(SOLIB)
else
 dox: $(DYLIB)
endif
	@echo
	@echo " Creating html documentation and logfile dohtml.log..."
	rm -f dohtml.log
	root -b -q dohtml.C 2>&1 >> dohtml.log | tee -a dohtml.log
	@echo " done."
	@echo

#clean:	rmcint rmobjs rmdep rmcore rmlib

mrproper:	$(MRPROPERS) rmbin rmbak rmbakmac rmhtml clean
	@echo " Done."
	@echo " "

tar:	mrproper
	@echo "Making tar-file"
	root -b -q -l -n tar.C
#	@tar cvf ../mars.tar --exclude=Root .rootrc *
#	@gzip -9 ../mars.tar

test:
	@echo "Calling make in mtest..."
	(cd mtest; $(MAKE)) 

testclean:
	@echo "Calling clean in mtest..."
	(cd mtest; $(MAKE) clean)
 
#Makefile.depend:
#	(! find ./ Makefile.depend -maxdepth 1 -empty 2> /dev/null && \
#	echo " Generating dependancies into Makefile.depend" && \
#	makedepend -- $(INCLUDES) -- $(PROGRAMS:=.cc) $(SRCS) $(SUBDIRS:=/*.cc) -w1024 -f- 2> /dev/null | grep -v Cint | grep -v "/usr/" > Makefile.depend && \
#	echo " ") || find -maxdepth 0 -true > /dev/null
#
#depend:	Makefile.depend	

# @endcode
