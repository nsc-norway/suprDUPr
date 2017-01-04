#include <utility>

#include <fstream>
#include <iostream>

#include <unordered_map>
#include <forward_list>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
//#include "gzip.hpp"
#include <boost/program_options.hpp>

#define BUFFER_SIZE 4*1024*1024
using namespace std;
namespace po = boost::program_options;

class InputSelector {
    // InputSelector class sets up the input stream based on 
    // command line arguments.
    private:
        istream* raw_input, *input_ptr;
        unique_ptr<istream> filtered_input;
        unique_ptr<char> buffer;
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
            
            char* buf_array = new char[BUFFER_SIZE];
            char* buf_array2 = new char[BUFFER_SIZE];
            raw_input->rdbuf()->pubsetbuf(buf_array, BUFFER_SIZE);
            buffer.reset(buf_array);

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
    unsigned long data[N] = {};
    
    TwoBitSequence(const char* seq, size_t str_len) {
        const unsigned long* blocks = (const unsigned long*) seq;
        for (size_t i=0; i<N; ++i) {
            data[i] = (blocks[i*4] & (0x0606060606060606ul)) >> 1;
            for (size_t j=0; j<3; ++j) {
                data[i] |= (blocks[i*4+j+1] & (0x0606060606060606ul)) << (j*2+1);
            }
        }
    }

    bool operator==(const TwoBitSequence& other) const {
        for (int i=0; i<N; ++i)
            if (data[i] != other.data[i])
                return false;
        return true;
    }
};

namespace std {
    template<size_t N>
    struct hash<TwoBitSequence<N>>
    {
        size_t operator()(const TwoBitSequence<N>& k) const
        {
            size_t hash = 0;
            for (int i=0; i<N; ++i) {
                hash ^= k.data[i];
            }
            return hash;
        }
    };
}


class Entry {
    public:
        int x, y;

        Entry(int x, int y) : x(x), y(y) {
        }
};

template<typename KEY>
class AnalysisHead {

    int winx, winy;
    size_t str_len;
    unsigned long duplicates = 0;

    public:
        typedef forward_list<Entry> entry_list;
        
        unordered_map<KEY, entry_list> data;

        AnalysisHead(size_t str_len, unsigned int winx, unsigned int winy) 
            : str_len(str_len), winx(winx), winy(winy) {
        }

        /* Enter a new sequence. Takes ownership of the string pointers.
         * This should be run inside an omp parallel block, as it creates tasks. */
        void enterPoint(int x, int y, const char* str) {
            entry_list& entries = data[KEY(str, str_len)];
            entries.emplace_front(x, y);
            entry_list::iterator it = entries.begin();
            ++it;
            for (; it != entries.end(); ++it) {
                Entry& entry = *it;
                if (abs(entry.x - x) < winx && abs(entry.y - y) < winy) {
                    duplicates++;
                }
            }
        }

        void endOfGroup() {
            data.clear();
        }

        unsigned long getTotal() {
            return duplicates;
        }
};

template <typename KEY>
pair<unsigned long, unsigned long> analysisLoop(
        int str_start, int str_len,
        int winx, int winy, istream& input,
        const string& header, const string& sequence
        ) {

    size_t start_to_coord_offset, start_to_y_coord_offset;
    size_t colons = 0, end_coords = 0; // temporary variables

    for (size_t i=0; i<header.size()-1; ++i) {
        if (header[i] == ':') {
            ++colons;
            if (colons == 5) start_to_coord_offset = i+1;
            if (colons == 6) start_to_y_coord_offset = i+1;
        }
        else if (header[i] == ' ') {
            end_coords = i;
        }
    }
    if (colons < 6 || !end_coords) {
        cerr << "Illumina format x/y coordinates not detected" << endl;
        return pair<unsigned long, unsigned long>(0, 1);
    }


    int x, y;
    x = atoi(header.c_str() + start_to_coord_offset);
    y = atoi(header.c_str() + start_to_y_coord_offset);
    size_t coord_to_seq_offset = header.size() - end_coords + 1;

    size_t i_record = 0, seq_len = sequence.size();

    // Read identifier prefix (before coordinates)
    char read_id[2][start_to_coord_offset];
    memcpy(read_id[0], &header[0], start_to_coord_offset);
    int ibuf = 1; // Identifier of where to write next read ID prefix (0/1)

    const size_t buf_size = ((str_len + 1) * 32 - 1) / 32; // Round up (pad zero)
    char* seq_data = new char[buf_size];
    sequence.copy(seq_data, str_len, str_start);

    cerr << "Started reading FASTQ file..." << endl;
    input.ignore(3 + seq_len); // Ignores "+\n" line & quality and \n

    AnalysisHead<KEY> analysisHead(seq_len, winx, winy);

    if (input) {
        analysisHead.enterPoint(x, y, seq_data);
        i_record++;
        while (input) { // Input loop
            input.read(read_id[ibuf], start_to_coord_offset);
            if (!input) break;
            if (memcmp(read_id[0], read_id[1], start_to_coord_offset) != 0) {
                analysisHead.endOfGroup();
            }
            ibuf ^= 1;

            input >> x;
            input.ignore(1);
            input >> y;
            input.ignore(coord_to_seq_offset + str_start);
            input.read(seq_data, str_len);
            size_t remain = seq_len - str_len - str_start;
            input.ignore(remain + 4 + seq_len); // Ignores "xxx\n+\nquality string\n"
            analysisHead.enterPoint(x, y, seq_data);
            if (++i_record % 1000000 == 0) cerr << "Read " << i_record << " records." << endl;
        }
    }
    return pair<unsigned long, unsigned long>(i_record, analysisHead.getTotal());
}



int main(int argc, char* argv[]) {

    string inputfile("-"), outputfile("-");
    unsigned int winx, winy;
    int first_base, last_base = -1;

    cerr << "Fastdup called with args: ";
    for (int i=1; i!=argc; ++i) {
        cerr << argv[i];
    }
    cerr << '\n';

    po::options_description visible("Allowed options");
    visible.add_options()
        ("winx,x", po::value<unsigned int>(&winx)->default_value(2500), "x coordinate window, +/- pixels")
        ("winy,y", po::value<unsigned int>(&winy)->default_value(2500), "y coordinate window, +/- pixels")
        ("start,s", po::value<int>(&first_base)->default_value(0),
            "First base position in reads to consider")
        ("end,e", po::value<int>(&last_base)->default_value(-1), 
            "Last base position in reads to consider (use -1 for all bases)")
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
    po::store(
            po::command_line_parser(argc, argv).options(all_options).positional(pos_desc).run(),
            vm
            );
    po::notify(vm);


    if (vm.count("help") > 0) {
        cerr << "usage: " << argv[0] << " [options] [input_file [output_file]] \n";
        cerr << visible << '\n';
        return 0;
    }
    
    InputSelector isel(inputfile);

    if (!isel.valid) {
        cerr << "Cannot open file " << inputfile << ": " << strerror(errno) << "\n";
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
    size_t str_start = min((size_t)first_base, seq_len);
    if (last_base == -1) last_base = seq_len;
    size_t str_len = min(seq_len - str_start, (size_t)last_base);

    cerr << "Using bases from " << str_start << " to " <<
        str_start+str_len << endl;

    pair<unsigned long, unsigned long> result;
    // Initial analysis to determine byte offsets (positions)
    if (str_len > 160) {
        result = analysisLoop<string>(
                str_start, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 128) {
        result = analysisLoop<TwoBitSequence<5>>(
                str_start, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 96) {
        result = analysisLoop<TwoBitSequence<4>>(
                str_start, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 64) {
        result = analysisLoop<TwoBitSequence<3>>(
                str_start, str_len, winx, winy, input, header, sequence
                );
    }
    else if (str_len > 32) {
        result = analysisLoop<TwoBitSequence<2>>(
                str_start, str_len, winx, winy, input, header, sequence
                );
    }
    else {
        result = analysisLoop<TwoBitSequence<1>>(
                str_start, str_len, winx, winy, input, header, sequence
                );
    }

    if (result.first == 0 && result.second != 0) {
        return 1; // error flag
    }
    else if (input.eof() && output.good()) {
        unsigned long records = result.first;
        unsigned long dupes = result.second;
        output <<  records << '\t' << dupes << endl;
        return 0;
    }
    else {
        if (input.bad()) {
            cerr << "Read error: " << strerror(errno) << endl;
        }
        else {
            cerr << "Unexpected problem!" << endl;
        }
        return 1;
    }
}



