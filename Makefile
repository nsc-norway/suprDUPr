all: suprDUPr suprDUPr.read_id

CFLAGS += -O3 -std=c++11

suprDUPr: suprDUPr.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -pthread -lz -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF)

suprDUPr.static: suprDUPr.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -static -pthread -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF) -lz

suprDUPr.read_id: suprDUPr.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -DOUTPUT_READ_ID -pthread -lz -lboost_program_options$(BOOST_LIB_SUFF) -lboost_iostreams$(BOOST_LIB_SUFF)

filterfq: filterfq.cpp
	$(CXX) -o $@ $^ $(CFLAGS) -pthread -lz -lboost_iostreams$(BOOST_LIB_SUFF)

duplicate-finder.subrange: duplicate-finder.subrange.cpp
	$(CXX) -Wall -o $@ $^ $(CFLAGS) -fopenmp -lz `ldconfig -p | awk -F' => ' '/ *libboost_iostreams\.so / { print $$2; }'`

clean:
	rm -f suprDUPr suprDUPr.read_id duplicate-finder.subrange
