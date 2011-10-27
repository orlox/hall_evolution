#if "t" (from toroidal) is defined, ignore the parts related to the poloidal part in the compilation
ifdef t
	TOROIDAL=-DTOROIDAL
endif
#if "o" (from ohm) is defined, ignore the parts related to the hall effect in the compilation
ifdef o
	PUREOHM=-DPUREOHM
endif
#if "s" (from simple) is defined, use s simple boundary condition, alpha=0, at the boundaries
ifdef s
	SIMPLE=-DSIMPLE
endif
#"all" compiles the code, ignoring the parts specified by TOROIDAL, PUREOHM and SIMPLE.
all:
	g++ -Wall -lgsl $(TOROIDAL) $(PUREOHM) $(SIMPLE) as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#"the same as all, but with debug symbols"
dall:
	g++ -g -Wall -lgsl as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
