all: fastdup fastdup.read_id

CPPFLAGS = -O3 -std=c++11

fastdup: fastdup.cpp
	$(CXX) -I /usr/include/boost148/ -L /usr/lib64/boost148/ -o $@ $^ $(CPPFLAGS) -lz -lboost_program_options -lboost_iostreams

fastdup.read_id: fastdup.read_id.cpp
	$(CXX) -o $@ $^ $(CPPFLAGS) -lz -lboost_program_options -lboost_iostreams

duplicate-finder.subrange: duplicate-finder.subrange.cpp
	$(CXX) -Wall -o $@ $^ $(CPPFLAGS) -fopenmp -lz `ldconfig -p | awk -F' => ' '/ *libboost_iostreams\.so / { print $$2; }'`

clean:
	rm -f fastdup fastdup.read_id duplicate-finder.subrange
