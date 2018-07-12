#include <fstream>
#include <iostream>
#include <iomanip>

#include <unordered_map>
#include <forward_list>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>

/*
 * filterdups.cpp
 *
 * Utility program to filter a FASTQ file based on ID strings from STDIN.
 *
 * Expected input is of the form outputted by suprDUPr.read_id:
 *
 * E00426:11:HCJFHALXX:1:1101:8724:7181	E00426:11:HCJFHALXX:1:1101:8329:6741
 * E00426:11:HCJFHALXX:1:1101:14022:7181	E00426:11:HCJFHALXX:1:1101:13768:7093
 * E00426:11:HCJFHALXX:1:1101:16579:7181	E00426:11:HCJFHALXX:1:1101:16589:7163
 *
 * The input may be tab-separated (as the suprDUPr.read_id output) or just contain one
 * column. FASTQ reads with ID matching strings from the first column are not written
 * to output. All other reads are written to the output. The input FASTQ file and the
 * read ID strings are expected to have the same sorting.
 *
 */

using namespace std;
namespace po = boost::program_options;

class InputSelector {
    // InputSelector class sets up the input stream from STDIN, or opens a file,
    // and detects whether the input is GZIP compressed.
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

int main(int argc, char* argv[]) {

    // Main function: Reads arguments and calls analysisLoop
    
    string inputfile("-"), outputfile("-");
    unsigned int winx, winy;
    int first_base, last_base = -1;
    size_t hash_bytes;
    bool region_sorted, unsorted;

    po::options_description visible("Allowed options");
    visible.add_options()
        ("winx,x", po::value<unsigned int>(&winx)->default_value(2500),
            "x coordinate window, +/- pixels")
        ("winy,y", po::value<unsigned int>(&winy)->default_value(2500),
            "y coordinate window, +/- pixels")
        ("start,s", po::value<int>(&first_base)->default_value(10),
            "First nucleotide position in reads to consider")
        ("end,e", po::value<int>(&last_base)->default_value(60), 
            "Last nucleotide position in reads to consider (use -1 for the end of the read)")
        ("region-sorted,r", po::bool_switch(&region_sorted),
            "Assume the input file is sorted by region (tile), but not by (y, x) coordinate "
            "within the region.")
        ("unsorted,u", po::bool_switch(&unsorted),
            "Process unsorted file (a large hash-size is recommended, see --hash-size). This "
            "mode requires all data to be stored in memory, and it is not well optimised.")
        ("hash-size", po::value<size_t>(&hash_bytes)->default_value(512*1024*8), 
            "Hash table size (bytes), must be a power of 2. (increase if winy>2500).")
        ("help,h", "Show this help message")
    ;
    po::options_description positionals("Positional options(hidden)");
    positionals.add_options()
        ("input-file", po::value<string>(&inputfile)->required(),
            "Input file, or - to read from STDIN")
        ("output-file", po::value<string>(&outputfile),
            "Output file, or - to write to STDOUT")
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
        if (vm.count("help") > 0) {
            printUsage(argv[0]);
            cerr << visible << '\n';
            cerr << "  Specify - for input_file to read from standard input.\n" << endl;
            return 0;
        }
        else if (vm.count("input-file") != 1) {
            cerr << "ERROR: Argument input_file is required." << endl;
            printUsage(argv[0]);
            cerr << visible << '\n';
            cerr << "  Specify - for input_file to read from stdin." << endl;
            return 1;
        }
        else {
            po::notify(vm);
        }
    }
    catch(po::error& e) 
    { 
      cerr << "ERROR: " << e.what() << "\n\n"; 
      printUsage(argv[0]);
      cerr << visible << endl; 
      return 1; 
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

    cerr << "-- suprDUPr v1.0 --\n";

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

    cerr << "Using nucleotide positions from " << first_base << " to "
         << first_base+str_len << endl;

    // Call the correct analysis loop for the specified string length
    // All these versions of the analysis loop are compiled as separate 
    // function, but only one is used for a given set of input parameters.
    Metrics result;
    if (str_len > 160) {
        cerr << "ERROR: Sorry, strings longer than 160 bytes are not supported "
            << "(use parameters --start, --end)" << endl;
        return 1;
    }
    else if (str_len > 128) {
        result = analysisLoop<TwoBitSequence<5>>(
                output,
                hash_bytes, first_base, str_len, winx, winy, region_sorted, unsorted,
                input, header, sequence
                );
    }
    else if (str_len > 96) {
        result = analysisLoop<TwoBitSequence<4>>(
                output,
                hash_bytes, first_base, str_len, winx, winy, region_sorted, unsorted,
                input, header, sequence
                );
    }
    else if (str_len > 64) {
        result = analysisLoop<TwoBitSequence<3>>(
                output,
                hash_bytes, first_base, str_len, winx, winy, region_sorted, unsorted,
                input, header, sequence
                );
    }
    else if (str_len > 32) {
        result = analysisLoop<TwoBitSequence<2>>(
                output,
                hash_bytes, first_base, str_len, winx, winy, region_sorted, unsorted,
                input, header, sequence
                );
    }
    else if (str_len > 0) {
        result = analysisLoop<TwoBitSequence<1>>(
                output,
                hash_bytes, first_base, str_len, winx, winy, region_sorted, unsorted,
                input, header, sequence
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
        cerr << "Completed. Analysed " << result.num_reads << " records." << endl;
#ifdef OUTPUT_READ_ID
        ostream& statsstream = cerr;
#else
        ostream& statsstream = output;
#endif
        statsstream << "NUM_READS\tREADS_WITH_DUP\tDUP_RATIO\n";
        statsstream << result.num_reads 
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
            cerr << "good=" << input.good() << ", eof=" << input.eof() << endl;
        }
        return 1;
    }
}



