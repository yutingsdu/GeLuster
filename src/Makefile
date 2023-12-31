#-- for now these MUST point to the included "samtools-0.x.x" and "gclib" sub-directories
BAM  := ./samtools-0.1.18
GDIR := ./gclib
#--

INCDIRS := -I. -I${GDIR} -I${BAM}

CXX   := $(if $(CXX),$(CXX),g++)

BASEFLAGS := -Wall -Wextra ${INCDIRS} -fsigned-char -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -std=c++0x -fno-strict-aliasing -fno-exceptions -fno-rtti
#for gcc 8+ add: -Wno-class-memaccess
GCCVER5 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 5)
ifeq "$(GCCVER5)" "1"
 BASEFLAGS += -Wno-implicit-fallthrough
endif

GCCVER8 := $(shell expr `${CXX} -dumpversion | cut -f1 -d.` \>= 8)
ifeq "$(GCCVER8)" "1"
  BASEFLAGS += -Wno-class-memaccess
endif

LINKER  := $(if $(LINKER),$(LINKER),g++)

LDFLAGS := $(if $(LDFLAGS),$(LDFLAGS),-g)

LDFLAGS += -L${BAM}

LIBS    := -lbam -lz

ifneq (,$(findstring nothreads,$(MAKECMDGOALS)))
 NOTHREADS=1 
endif

#detect MinGW (Windows environment)
ifneq (,$(findstring mingw,$(shell ${CXX} -dumpmachine)))
 WINDOWS=1
endif

# Misc. system commands
ifdef WINDOWS
 RM = del /Q
else
 RM = rm -f
endif

# File endings
ifdef WINDOWS
 EXE = .exe
else
 EXE =
endif

# Non-windows systems need pthread
ifndef WINDOWS
 ifndef NOTHREADS
   LIBS := -pthread ${LIBS}
   BASEFLAGS += -pthread
 endif
endif

ifdef NOTHREADS
  BASEFLAGS += -DNOTHREADS
endif

DMACH := $(shell ${CXX} -dumpmachine)

ifneq (,$(filter %release %static, $(MAKECMDGOALS)))
  # -- release build
  RELEASE_BUILD=1
  CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O3)
  CXXFLAGS += -DNDEBUG $(BASEFLAGS)
else
  ifneq (,$(filter %memcheck %memdebug %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
     #use sanitizer in gcc 4.9+
     GCCVER49 := $(shell expr `${CXX} -dumpversion | cut -f1,2 -d.` \>= 4.9)
     ifeq "$(GCCVER49)" "0"
       $(error gcc version 4.9 or greater is required for this build target)
     endif
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     SANLIBS :=
     ifneq (,$(filter %tsan %tcheck %thrcheck, $(MAKECMDGOALS)))
        # thread sanitizer only (incompatible with address sanitizer)
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=thread -fsanitize=undefined $(BASEFLAGS)
        SANLIBS := -ltsan
     else
        # address sanitizer
        CXXFLAGS += -fno-omit-frame-pointer -fsanitize=undefined -fsanitize=address $(BASEFLAGS)
        SANLIBS := -lasan
     endif
     ifeq "$(GCCVER5)" "1"
       CXXFLAGS += -fsanitize=bounds -fsanitize=float-divide-by-zero -fsanitize=vptr
       CXXFLAGS += -fsanitize=float-cast-overflow -fsanitize=object-size
       #CXXFLAGS += -fcheck-pointer-bounds -mmpx
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG -fno-common -fstack-protector
     LIBS := ${SANLIBS} -lubsan -ldl ${LIBS}
  else
     #just plain debug build
     DEBUG_BUILD=1
     CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-g -O0)
     #CXXFLAGS := $(if $(CXXFLAGS),$(CXXFLAGS),-ggdb -g3 -O0 -fvar-tracking-assignments -fno-omit-frame-pointer)
     ifneq (, $(findstring darwin, $(DMACH)))
        CXXFLAGS += -gdwarf-3
     endif
     CXXFLAGS += -DDEBUG -D_DEBUG -DGDEBUG $(BASEFLAGS)
  endif
endif

ifdef RELEASE_BUILD
 ifneq (,$(findstring static, $(MAKECMDGOALS)))
    STATIC_CLIB=1
 endif
endif

ifdef STATIC_CLIB
 LDFLAGS += -static-libgcc -static-libstdc++
endif

ifdef DEBUG_BUILD
  #$(warning Building DEBUG version of get_alignment_info.. )
  DBG_WARN=@echo
  DBG_WARN+='WARNING: built DEBUG version [much slower], use "make clean release" for a faster, optimized version of the program.'
endif

OBJS := ${GDIR}/GBase.o ${GDIR}/GArgs.o ${GDIR}/GStr.o ${GDIR}/GBam.o \
 ${GDIR}/gdna.o ${GDIR}/codons.o ${GDIR}/GFastaIndex.o ${GDIR}/GFaSeqGet.o ${GDIR}/gff.o 

ifneq (,$(filter %memtrace %memusage %memuse, $(MAKECMDGOALS)))
    CXXFLAGS += -DGMEMTRACE
    OBJS += ${GDIR}/proc_mem.o
endif

ifndef NOTHREADS
 OBJS += ${GDIR}/GThreads.o 
endif

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

#OBJS += rlink.o tablemaker.o tmerge.o common.o
OBJS += tmerge.o

all release static debug: get_alignment_info${EXE} PseudoRef${EXE} get_cluster${EXE} FinalCluster${EXE} MergeFiles${EXE}
memcheck memdebug tsan tcheck thrcheck: get_alignment_info${EXE} PseudoRef${EXE} get_cluster${EXE} FinalCluster${EXE} MergeFiles${EXE}
memuse memusage memtrace: get_alignment_info${EXE} PseudoRef${EXE} get_cluster${EXE} FinalCluster${EXE} MergeFiles${EXE}
nothreads: get_alignment_info${EXE} PseudoRef${EXE} get_cluster${EXE} FinalCluster${EXE} MergeFiles${EXE}

${GDIR}/GBam.o : $(GDIR)/GBam.h
get_alignment_info.o : $(GDIR)/GBitVec.h $(GDIR)/GHash.hh $(GDIR)/GBam.h
#rlink.o : rlink.h tablemaker.h $(GDIR)/GBam.h $(GDIR)/GBitVec.h 
tmerge.o :  tmerge.h
#tablemaker.o : tablemaker.h rlink.h
#common.o : common.h
utility.o : utility.h
${BAM}/libbam.a: 
	cd ${BAM} && make lib
get_alignment_info: ${BAM}/libbam.a $(OBJS) get_alignment_info.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
PseudoRef: PseudoRef.o utility.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
get_cluster: get_cluster.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
FinalCluster: FinalCluster.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
MergeFiles: MergeFiles.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LIBS}
	@echo
	${DBG_WARN}
test demo tests: get_alignment_info${EXE}
	@./run_tests.sh
.PHONY : clean cleanall cleanAll allclean

# target for removing all object files

#	echo $(PATH)
clean:
	${RM} get_alignment_info${EXE} get_alignment_info.o* PseudoRef${EXE} PseudoRef.o* get_cluster${EXE} get_cluster.o*  FinalCluster${EXE} FinalCluster.o* MergeFiles${EXE} MergeFiles.o* $(OBJS)
	${RM} core.*
	${RM} utility.o FinalCluster.o 
allclean cleanAll cleanall:
	cd ${BAM} && make clean
	${RM} get_alignment_info${EXE} get_alignment_info.o* PseudoRef${EXE} PseudoRef.o* get_cluster${EXE} get_cluster.o*  FinalCluster${EXE} FinalCluster.o* MergeFiles${EXE} MergeFiles.o* $(OBJS) 
	${RM} core.*
	${RM} utility.o
