all: suprDUPr suprDUPr.read_id filterfq

CFLAGS += -O3 -std=c++11

suprDUPr: suprDUPr.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -pthread -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz

suprDUPr.read_id: suprDUPr.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -DOUTPUT_READ_ID -pthread -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz

filterfq: filterfq.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -pthread -lboost_iostreams$(BOOST_LIB_SUFF) -lz

duplicate-finder.subrange: duplicate-finder.subrange.cpp
	$(CXX) -Wall -o $@ $^ $(CFLAGS) -fopenmp -lz `ldconfig -p | awk -F' => ' '/ *libboost_iostreams\.so / { print $$2; }'`

clean:
	rm -f suprDUPr suprDUPr.read_id duplicate-finder.subrange filterfq
