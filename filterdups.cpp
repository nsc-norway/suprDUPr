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


int main(int argc, char* argv[]) {

    // Main function: Reads arguments and calls analysisLoop
    string inputfile("-");

    istream* raw_input, *input_ptr;
    unique_ptr<istream> filtered_input;
    boost::iostreams::filtering_istream in;

    bool valid;

    // Disable sync with printf, etc.
    ios_base::sync_with_stdio(false);
    // Use a simple locale, for speed
    setlocale(LC_ALL,"C");
    // Prevents flushing cout when reading from cin
    cin.tie(nullptr);

    istream* inputfile;
    ifstream file_input;
    file_input.open(filename, ios_base::in | ios_base::binary);
    if (file_input.fail()) {
        valid = false;
        return;
    }
    raw_input = &file_input;
        
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

    cerr << "-- filterdups v1.0 --\n";

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


    return 0;
}



