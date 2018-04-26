all: fastdup fastdup.read_id

LOCALFLAGS = -O3 -std=c++11

fastdup: fastdup.cpp
	$(CXX) -o $@ $^ $(LOCALFLAGS) -lz -lboost_program_options -lboost_iostreams

fastdup.read_id: fastdup.cpp
	$(CXX) -o $@ $^ $(LOCALFLAGS) -DOUTPUT_READ_ID -lz -lboost_program_options -lboost_iostreams

duplicate-finder.subrange: duplicate-finder.subrange.cpp
	$(CXX) -Wall -o $@ $^ $(LOCALFLAGS) -fopenmp -lz `ldconfig -p | awk -F' => ' '/ *libboost_iostreams\.so / { print $$2; }'`

clean:
	rm -f fastdup fastdup.read_id duplicate-finder.subrange
