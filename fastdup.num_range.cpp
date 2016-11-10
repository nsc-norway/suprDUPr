#include <utility>

/**
 * Getting the number of entries "in range", i.e. for an y-window of 2500,
 * how many reads are located in the y-window at any time.
 *
 * Result: approx. 170k entries
 *
 * The purpose of this calculation is to figure out what kind of data 
 * structure may be used. E.g. a 400k bucket hash table would be large
 * compared to CPU cache, which could be a problem.
 */

//#include <unordered_map>
#include <fstream>
#include <iostream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
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
                }
                raw_input = &file_input;
            }
            
            char* buf_array = new char[BUFFER_SIZE];
            raw_input->rdbuf()->pubsetbuf(buf_array, BUFFER_SIZE);
            buffer.reset(buf_array);

            uint16_t firstInRow_two;
            uint8_t byte1, byte2;
            (*raw_input) >> byte1 >> byte2;
            if ( !(*raw_input).good() ) {
                valid = false;
                return;
            }
            raw_input->putback(byte2);
            raw_input->putback(byte1);

            if (byte1 == 0x1f && byte2 == 0x8b) {
                in.push(boost::iostreams::gzip_decompressor());
                in.push(*raw_input);
                input = &in;
            }
            else {
                input = raw_input;
            }
        }
};


class Entry {
    private:
        Entry(Entry& source) {}
    public:
        int x, y;
        unsigned int duplicates = 0;

        // Entry object takes ownership of the string pointers and
        // releases them when destroyed
        Entry() : x(-1), y(-1) {
        }

        Entry(Entry&& source) noexcept
           : x(source.x), y(source.y), duplicates(source.duplicates) {
        }

};


class DnaHasher {
    
    public:
        size_t len;

        DnaHasher(size_t len) : len(len) {
        }

        long hash(const char* seq) {
            long val = 0;
            //#pragma omp simd
            for (int i=0; i<len; ++i) {
                val |= ((seq[i] >> 1) & 3) << (i*2);
            }
            return val;
        }
};


class AnalysisHead {

    int winx, winy;
    size_t str_len;
    unsigned long duplicates = 0;

    public:
        
        hash<string> str_hash;
        //unordered_map<string, Entry> data;
        DnaHasher hasher;
        list<pair<int, int>> range;
        int n_in_range = 0, n_in_row = 0, current_y = 0;

        AnalysisHead(size_t str_len, unsigned int winx, unsigned int winy) 
            : str_len(str_len), winx(winx), winy(winy), hasher(str_len) {
        }

        /* Enter a new sequence. Takes ownership of the string pointers.
         * This should be run inside an omp parallel block, as it creates tasks. */
        void enterPoint(int x, int y, const char* seq) {
            //long test = hasher.hash(seq);
            if (y != current_y) {
                range.emplace_back(current_y, n_in_row);
                n_in_row = 0;
                while (!range.empty()) {
                    if (range.front().first < y - winy) {
                        n_in_range -= range.front().second;
                        range.pop_front();
                    }
                    else {
                        break;
                    }
                }
                current_y = y;
                cout << n_in_range << '\n';
            }
            n_in_row++;
            n_in_range++;
            //size_t h = str_hash(string(seq, str_len));
            //Entry& entry = data[string(seq, str_len)];
            //if (entry.x != -1) {
            //    duplicates++;
            //}
            //else {
            //}
        }

        void endOfGroup() {
            //data.clear();
            range.clear();
            n_in_range = 0;
            n_in_row = 0;
        }

        pair<unsigned int, unsigned int> getTotal() {
            return pair<unsigned int, unsigned int>(0, 0);
        }

};



int main(int argc, char* argv[]) {

    string inputfile("-"), outputfile("-");
    unsigned int winx, winy;
    int first_base, last_base = -1;
    po::options_description visible("Allowed options");
    visible.add_options()
        ("winx,x", po::value<unsigned int>(&winx)->default_value(2500), "x coordinate window, +/- pixels")
        ("winy,y", po::value<unsigned int>(&winy)->default_value(2500), "y coordinate window, +/- pixels")
        ("start,s", po::value<int>(&first_base)->default_value(25),
            "First base position in reads to consider")
        ("end,e", po::value<int>(&last_base)->default_value(75), 
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
    int i;
    size_t start_to_coord_offset = 0, coord_to_seq_offset = 0, seq_str_len = 0;
    string header, sequence;
    getline(input, header);
    getline(input, sequence);
    int colons = 0, end_coords = 0, start_to_y_coord_offset;
    for (int i=0; i<header.size()-1; ++i) {
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
        return 1;
    }
    int x, y;
    x = atoi(header.c_str() + start_to_coord_offset);
    y = atoi(header.c_str() + start_to_y_coord_offset);
    coord_to_seq_offset = header.size() - end_coords + 1;

    size_t i_record = 0, seq_len = sequence.size();
    // Read index within the sequence string (i.e. nucleotide position)
    size_t str_start = min((size_t)first_base, seq_len);
    if (last_base == -1) last_base = seq_len;
    size_t str_len = min(seq_len - str_start, (size_t)last_base);

    // Read identifier prefix (before coordinates)
    char read_id[2][start_to_coord_offset];
    memcpy(read_id[0], &header[0], start_to_coord_offset);
    int ibuf = 1; // Identifier of where to write next read ID prefix (0/1)

    char* seq_data = new char[str_len];
    sequence.copy(seq_data, str_len, str_start);

    cerr << "Started reading FASTQ file..." << endl;
    input.ignore(3 + seq_len); // Ignores "+\n" line & quality and \n

    AnalysisHead analysisHead(seq_len, winx, winy);
    if (input) {
        analysisHead.enterPoint(x, y, seq_data);
        i_record++;
        while (true) {
            input.read(read_id[ibuf], start_to_coord_offset);
            if (!input) break;
            if (memcmp(read_id[0], read_id[1], start_to_coord_offset) != 0) {
                cerr << "END OF GROUP!!" << endl;
                cerr << "\"";
                cerr.write(read_id[0], start_to_coord_offset);
                cerr << "\" vs \"";
                cerr.write(read_id[1], start_to_coord_offset);
                cerr << "\"\n";
                analysisHead.endOfGroup();
            }
            ibuf ^= 1;

            // Input loop
            input >> x;
            input.ignore(1);
            input >> y;
            //seq_data = new char[str_len];

            input.ignore(coord_to_seq_offset + str_start);
            input.read(seq_data, str_len);
            size_t remain = seq_len - str_len - str_start;
            input.ignore(remain + 4 + seq_len); // Ignores "xxx\n+\nquality string\n"
            analysisHead.enterPoint(x, y, seq_data);
            if (++i_record % 1000000 == 0) cerr << "Read " << i_record << " records." << endl;
        }
    }
    pair<unsigned int, unsigned int> totals = analysisHead.getTotal();
    cout << "Number of reads " << i_record << endl;
    cout << "Number of duplicates " << totals.first << endl;
    cout << "Number of duplicates (external) " << totals.second << endl;

    return 0;
}

