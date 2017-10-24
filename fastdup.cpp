#include <utility>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <unordered_map>
#include <forward_list>

#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
#include "gzip.hpp"
#include <boost/program_options.hpp>

#define MAX_LEN 256

using namespace std;
namespace po = boost::program_options;

class InputSelector {
    // InputSelector class sets up the input stream to STDIN or opens a file,
    // and detects GZIP compressed data.
    private:
        istream* raw_input, *input_ptr;
        unique_ptr<istream> filtered_input;
        ifstream file_input;
        boost::iostreams::filtering_istream in;

    public:
        istream* input;
        bool valid;

        InputSelector(const string& filename) {
            // Disable sync with printf, etc.
            ios_base::sync_with_stdio(false);
            // Use a simple locale, for speed
            setlocale(LC_ALL,"C");
            if (filename == "-") {
                raw_input = &cin;
                // Prevents flushing cout when reading from cin
                cin.tie(nullptr);
            }
            else {
                file_input.open(filename, ios_base::in | ios_base::binary);
                if (file_input.fail()) {
                    valid = false;
                    return;
                }
                raw_input = &file_input;
            }
            
            uint8_t byte1, byte2;

            (*raw_input) >> byte1 >> byte2;

            if ( !raw_input->good() ) {
                valid = false;
                return;
            }
            raw_input->putback(byte2);
            raw_input->putback(byte1);

            if (byte1 == 0x1f && byte2 == 0x8b) {
                cerr << "Input is GZIP compressed" << endl;
                in.push(boost::iostreams::gzip_decompressor());
                in.push(*raw_input);
                input = &in;
            }
            else {
                input = raw_input;
            }
            valid = true;
        }
};

template <size_t N>
class TwoBitSequence {
public:
    array<unsigned long, N> data = {};
    // Have half the number of unk's for N bases, but round up
    array<unsigned long, (N+1)/2> unk = {};
    
    inline TwoBitSequence(const char* seq) {
        /* Conversion of ASCII string to 2-bit codes + unknown (N) flag:
         *        vv
         *  A 1000001
         *  C 1000011
         *  G 1000111
         *  T 1010100
         *  N 1001110
         *
         * So we use a bit mask of 0x06 = 00000110 to fish out these bits, and
         * then shift them as appropriate.
         *
         * The 2-bit codes aren't stored in order, making it a bit more efficient
         * to read them in 8-byte chunks. One can fit 4 2-bit bases into the space
         * of a single ASCII character (byte). The bases are stored in an unique
         * order (order is irrelevant as long as equality is preserved):
         * 1,9,17,25,2,10,18,26,3,11,19,27,.....,8,16,24,32
         *
         * The sequence length is not needed. It is the responsibility of the 
         * caller to call the TwoBitSequence<N> with the correct N, where N
         * is seq_len/32 rounded up. Also, the buffer passed to the constructor
         * must be at least N*32 bytes long, zero padded if necessary.
         *
         * The handling of N was added later.
         */
        const unsigned long* blocks = (const unsigned long*) seq;
        for (size_t i=0; i<N; ++i) {
            data[i] = (blocks[i*4] & (0x0606060606060606ul)) >> 1;
            size_t j;
            for (j=0; j<3; ++j) {
                data[i] |= (blocks[i*4+j+1] & (0x0606060606060606ul)) << (j*2+1);
            }
            size_t jshift, nshift;
	    j = 0;
            if ((i & 1) == 1) 
                nshift = 4;
            else
                nshift = 0;
	    jshift = nshift;
            for (; jshift<3; ++jshift) {
                unk[i/2] |= (blocks[i*4+(j++)] & (0x0808080808080808ul)) >> (3-jshift);
            }
            for (; jshift<4+nshift; ++jshift) {
                unk[i/2] |= (blocks[i*4+(j++)] & (0x0808080808080808ul)) << (jshift-3);
            }
        }
    }

    inline bool operator==(const TwoBitSequence& other) const {
        return data == other.data && unk == other.unk;
    }

    inline size_t hash() const {
        // Hash ignores unk; treats N as G. As N is uncommon, it's 
        // not worth the effort.
        size_t hash = 0;
        for (size_t i=0; i<N; ++i) {
            hash += data[i];
        }
        return hash;
    }
};

template<typename VALUE>
class Entry {
    // This class represents a single sequence read at a specific position 
    // on the flow cell. 
    public:
        Entry* next = nullptr;
        short group;
        int x, y;
        VALUE value;

        Entry(short group, int x, int y, VALUE value) :
            group(group), x(x), y(y), value(value) {
        }

        Entry(short group, int x, int y, const char* seq) :
            group(group), x(x), y(y), value(seq) {
        }
};

class Metrics {
    // Collects totals
    public:
        bool error = false;
        unsigned long reads_with_duplicates = 0;
        unsigned long num_reads = 0;
};

template<typename VALUE>
class AnalysisHead {

    /* The AnalysisHead class receives read one by one from the main program, and 
     * manages the processing of rows, and groups (tiles). */

    const size_t hash_size;
    const size_t mask;
    int winx, winy;
    size_t str_len;

    short group = 0;

    typedef Entry<VALUE> Ent;
    Ent** data = nullptr;

    public:
    Metrics metrics;
        
        AnalysisHead(size_t hash_bytes, unsigned int winx, unsigned int winy) 
            : hash_size(hash_bytes/sizeof(Ent*)), mask(hash_size-1),
                winx(winx), winy(winy) {
            data = new Ent*[hash_size];
        }

        ~AnalysisHead() {
            delete data;
        }

        void enterPoint(int x, int y, const char* str) {
            Ent* new_entry = new Ent(group,x,y,str);
            Ent** entry_ptr = &data[new_entry->value.hash() & mask];
            unsigned int any_duplicate_found = 0;
            while (*entry_ptr) {
                Ent* entry = (*entry_ptr);
                if ((y - entry->y) > winy || entry->group != group) {
                    *entry_ptr = entry->next;
                    delete entry;
                }
                else {
                    if (any_duplicate_found == 0
                            && abs(entry->x - x) < winx
                            && entry->value == new_entry->value) {
                        any_duplicate_found = 1;
                        (new_entry)->next = entry;
                        break;
                    }
                    entry_ptr = &entry->next;
                }
            }
            metrics.reads_with_duplicates += any_duplicate_found;
            metrics.num_reads++;
            *entry_ptr = new_entry;
        }

        void endOfGroup() {
            group++;
        }
};

template <typename VALUE>
Metrics analysisLoop(
        size_t hash_bytes, int str_start, int str_len,
        int winx, int winy, istream& input,
        const string& header, const string& sequence
        ) {

    /*
     * Function analysisLoop is called by main program to run the actual analysis.
     *
     * This is separated into a different function in order to be able to use a 
     * type parameter (VALUE), so it can call the corresponding AnalysisHead efficiently.
     */

    // Do some more work to determine the header format. In Illumina format, the 
    // x coordinate is after the fifth colon and the y coordinate after the sixth 
    // colon up to a space.
    size_t start_to_coord_offset = 0, start_to_y_coord_offset = 0;
    size_t colons = 0, end_coords = 0; // temporary variables
    size_t i;
    for (i=0; i<header.size(); ++i) {
        if (header[i] == ':') {
            ++colons;
            if (colons == 5) start_to_coord_offset = i+1;
            if (colons == 6) start_to_y_coord_offset = i+1;
        }
        else if (header[i] == ' ') {
            end_coords = i;
        }
    }
    if (end_coords == 0)
        end_coords = i;
    if (colons < 6) {
        cerr << "Illumina format x/y coordinates not detected" << endl;
        Metrics m;
        m.error = true;
        return m;
    }

    // Read the x/y of the first entry in a special way, because we have already
    // read the lines into memory.
    int x, y;
    x = atoi(header.c_str() + start_to_coord_offset);
    y = atoi(header.c_str() + start_to_y_coord_offset);
    size_t coord_to_seq_offset = header.size() - end_coords + 1;

    size_t seq_len = sequence.size();
    if (seq_len >= MAX_LEN) {
        cerr << "Sequence is too long, max supported is: " << MAX_LEN << "." << endl;
        Metrics m;
        m.error = true;
        return m;
    }

    // Read identifier prefix (before coordinates)
    char read_id[2][start_to_coord_offset];
    memcpy(read_id[0], &header[0], start_to_coord_offset);
    int ibuf = 1; // Identifier of where to write next read ID prefix (0/1)

    const size_t pad_to = sizeof(unsigned long) * 4;
    const size_t buf_size = ((str_len + 1) * pad_to - 1) / pad_to; // Round up (pad zero)
    char* seq_data = new char[buf_size];
    char* linebuf = new char[MAX_LEN];
    char* dummybuf = new char[MAX_LEN];
    sequence.copy(seq_data, str_len, str_start);

    // We also accept a header between sequence and quality, doesn't need to be just
    // "+". For example, some files from SRA have this.
    string middle_header;
    getline(input, middle_header);

    if (middle_header.size() == 0 || middle_header[0] != '+') {
        cerr << "The line after the sequence doesn't conform to the expected " 
            << "format." << endl;
        Metrics m;
        m.error = true;
        return m;
    }

    input.ignore(1 + seq_len); // Ignores quality scores and \n

    AnalysisHead<VALUE> analysisHead(hash_bytes, winx, winy);

    cerr << "Started reading FASTQ file..." << endl;

    if (input) {
        // This point was already read
        analysisHead.enterPoint(x, y, seq_data);

        while (input) { // Input loop
            // Read the read ID (query name) up to the coordinates
            input.read(read_id[ibuf], start_to_coord_offset);
            if (!input) break;
            // If it doesn't match the last one, signal "end of group" (tile)
            if (memcmp(read_id[0], read_id[1], start_to_coord_offset) != 0) {
                analysisHead.endOfGroup();
            }
            if (read_id[ibuf][0] != '@') {
                cerr << "Header format error: expected it to start with '@', got '" 
                     << read_id[ibuf][0] << "' at record #" << analysisHead.metrics.num_reads
                     << "." << endl;
                Metrics m;
                m.error = true;
                return m;
            }
            ibuf = ibuf == 0 ? 1 : 0;

            // Read the coordinates, then ignore the rest of the header line
            char colon_test;
            input >> x >> colon_test >> y;
            input.ignore(coord_to_seq_offset);

            if (colon_test != ':') {
                cerr << "Invalid file format detected. All reads must be of the same length, "
                     << "and the header must be the standard Illumina header." << endl;
                cerr << "Expected ':' at line " << analysisHead.metrics.num_reads * 4 + 1
                     << ", got '" << (colon_test ? colon_test : '0') << "' (with x=" << x
                     << ", y=" << y << ")." << endl;
                Metrics m;
                m.error = true;
                return m;
            }

            // Read sequence substring
            input.getline(linebuf, MAX_LEN);

            // Number of characters read including end of line
            size_t num_read = input.gcount();
            // Ignore quality header, then the quality, to the end of the record
            // Quality header length is unspecified, need to use getline
            input.getline(dummybuf, MAX_LEN);
            input.ignore(num_read);

            if (num_read >= 1 + str_len + str_start) {
                memcpy(seq_data, linebuf + str_start, str_len);
                analysisHead.enterPoint(x, y, seq_data);
            }

            if (analysisHead.metrics.num_reads % 1000000 == 0)
                cerr << "Analysed " << setw(9) << analysisHead.metrics.num_reads
                    << " reads." << endl;
        }
    }
    return analysisHead.metrics;
}



int main(int argc, char* argv[]) {

    // Main function: Reads arguments and calls analysisLoop
    
    string inputfile("-"), outputfile("-");
    unsigned int winx, winy;
    int first_base, last_base = -1;
    size_t hash_bytes;

    po::options_description visible("Allowed options");
    visible.add_options()
        ("winx,x", po::value<unsigned int>(&winx)->default_value(2500), "x coordinate window, +/- pixels")
        ("winy,y", po::value<unsigned int>(&winy)->default_value(2500), "y coordinate window, +/- pixels")
        ("start,s", po::value<int>(&first_base)->default_value(10),
            "First base position in reads to consider")
        ("end,e", po::value<int>(&last_base)->default_value(60), 
            "Last base position in reads to consider (use -1 for all bases)")
        ("hash-size", po::value<size_t>(&hash_bytes)->default_value(512*1024*8), 
            "Hash table size (bytes), *must* be a power of 2. (increase if winy>2500).")
        ("help,h", "Show this help message")
    ;
    po::options_description positionals("Positional options(hidden)");
    positionals.add_options()
        ("input-file", po::value<string>(&inputfile), "Input file, or - to read from STDIN")
        ("output-file", po::value<string>(&outputfile), "Output file, or - to write to STDOUT")
    ;
    po::options_description all_options("Allowed options");
    all_options.add(visible);
    all_options.add(positionals);

    po::positional_options_description pos_desc;
    pos_desc.add("input-file", 1);
    pos_desc.add("output-file", 1);

    po::variables_map vm;
    try {
        po::store(
                po::command_line_parser(argc, argv).options(all_options).positional(pos_desc).run(),
                vm
                );
        po::notify(vm);
    }
    catch(po::error& e) 
    { 
      cerr << "ERROR: " << e.what() << "\n\n"; 
      cerr << visible << endl; 
      return 1; 
    } 


    if (vm.count("help") > 0) {
        cerr << "usage: " << argv[0] << " [options] [input_file [output_file]] \n";
        cerr << visible << '\n';
        return 0;
    }
    
    InputSelector isel(inputfile);

    if (!isel.valid) {
        if (inputfile == "-") {
            cerr << "ERROR: Cannot open standard input: " << strerror(errno) << "\n";
        }
        else {
            cerr << "ERROR: Cannot open file " << inputfile << ": " << strerror(errno) << "\n";
        }
        return 1;
    }

    ostream* output_ptr = &cout;
    ofstream output_file;
    if (vm.count("output-file") > 0) {
        output_file.open(vm["output-file"].as<string>(), ios_base::out);
        output_ptr = &output_file;
    }

    istream&input = *isel.input;
    ostream&output = *output_ptr;

    string header, sequence;

    getline(input, header);
    getline(input, sequence);

    size_t seq_len = sequence.size();
    // Read index within the sequence string (i.e. nucleotide position)
    if (last_base == -1) last_base = seq_len;
    size_t str_len = min(seq_len - first_base, (size_t)(last_base - first_base));

    cerr << "Using bases from " << first_base << " to " << first_base+str_len << endl;

    // Call the correct analysis loop for the specified string length
    Metrics result;
    if (str_len > 160) {
        cerr << "ERROR: Sorry, strings longer than 160 bytes are not supported "
            << "(use parameters --start, --end)" << endl;
        return 1;
    }
    else if (str_len > 128) {
        result = analysisLoop<TwoBitSequence<5>>(
                hash_bytes, first_base, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 96) {
        result = analysisLoop<TwoBitSequence<4>>(
                hash_bytes, first_base, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 64) {
        result = analysisLoop<TwoBitSequence<3>>(
                hash_bytes, first_base, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 32) {
        result = analysisLoop<TwoBitSequence<2>>(
                hash_bytes, first_base, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 0) {
        result = analysisLoop<TwoBitSequence<1>>(
                hash_bytes, first_base, str_len, winx, winy, input, header, sequence
                );
    }
    else {
        cerr << "ERROR: Zero length sequence to compare. Perhaps the format of "
            << "the file was not understood."<< endl;
        return 1;
    }

    if (result.error) {
        return 1; // error flag
    }
    else if (input.eof() && output.good()) {
	cerr << "Completed. Processed " << result.num_reads << " records." << endl;
        output << "NUM_READS\tREADS_WITH_DUP\tDUP_RATIO\n";
        output << result.num_reads 
            << '\t' << result.reads_with_duplicates 
            << '\t' << result.reads_with_duplicates * 1.0 / result.num_reads
            << endl;
        return 0;
    }
    else {
        if (input.bad()) {
            cerr << "ERROR: read: " << strerror(errno) << endl;
        }
        else {
            cerr << "ERROR: Unexpected problem!" << endl;
        }
        return 1;
    }
}



