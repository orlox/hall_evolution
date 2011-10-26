#"all" compiles the code including poloidal and toroidal fields and hall terms
all:
	g++ -Wall -lgsl as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#"t" compiles the code using only the toroidal field, and includes hall terms
t:
	g++ -Wall -lgsl -DTOROIDAL as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#"o" compiles the code including the poloidal and toroidal fields, but without hall terms (only resistivity)
o:
	g++ -Wall -lgsl -DPUREOHM as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#"to" compiles the code using only the toroidal field, and ignoring the hall terms (only resistivity)
to:
	g++ -Wall -lgsl -DTOROIDAL -DPUREOHM as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
#Rules that start with "d" are equivalent to the previous one, but include debugging symbols
dall:
	g++ -g -Wall -lgsl as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
dt:
	g++ -g -Wall -lgsl -DTOROIDAL as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
do:
	g++ -g -Wall -DPUREOHM as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
dto:
	g++ -g -Wall -DTOROIDAL -DPUREOHM as_full.cpp io.cpp sim.cpp initial.cpp -o as_full_run
