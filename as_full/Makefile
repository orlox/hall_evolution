#if "t" (from toroidal) is defined, ignore the parts related to the poloidal part in the compilation
ifdef t
	TOROIDAL=-DTOROIDAL
endif
#if "o" (from ohm) is defined, ignore the parts related to the hall effect in the compilation
ifdef o
	PUREOHM=-DPUREOHM
endif
#if mp is defined, use openmp
ifdef mp
	OPENMP=-fopenmp
endif
#if sc is defined, use supercondunctor BC
ifdef sc
	SC=-DSC
endif
#"all" compiles the code without debugging symbols.
all:
	g++ $(OPENMP) -Wall -lgsl $(TOROIDAL) $(PUREOHM) $(SC) as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#"the same as all, but with debug symbols"
dall:
	g++ $(OPENMP) -g -Wall -lgsl as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
