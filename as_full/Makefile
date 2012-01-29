#if "t" (from toroidal) is defined, ignore the parts related to the poloidal part in the compilation
ifdef t
	TOROIDAL=-DTOROIDAL
endif
#if "o" (from ohm) is defined, ignore the parts related to the hall effect in the compilation
ifdef o
	PUREOHM=-DPUREOHM
endif
#"all" compiles the code without debugging symbols.
all:
	g++ -Wall -lgsl $(TOROIDAL) $(PUREOHM) as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#"the same as all, but with debug symbols"
dall:
	g++ -g -Wall -lgsl as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
