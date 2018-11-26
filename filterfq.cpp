#include <fstream>
#include <iostream>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

/*
 * filterdups.cpp
 *
 * Utility program to filter a FASTQ file based on ID strings from STDIN. The strings
 * should not include a leading '@' character.
 *
 *  STDIN: Read-IDs --\
 *                     > filterdups -> STDOUT: filtered FASTQ data
 *  Input FASTQ file _/
 *
 *
 * Expected STDIN input is of the form outputted by suprDUPr.read_id:
 *
 * E00426:11:HCJFHALXX:1:1101:8724:7181	E00426:11:HCJFHALXX:1:1101:8329:6741
 * E00426:11:HCJFHALXX:1:1101:14022:7181	E00426:11:HCJFHALXX:1:1101:13768:7093
 * E00426:11:HCJFHALXX:1:1101:16579:7181	E00426:11:HCJFHALXX:1:1101:16589:7163
 *
 * The input may be tab-separated (as the suprDUPr.read_id output) or just contain one
 * column. FASTQ reads with ID matching strings from the first column are not written
 * to output. All other reads are written to the output. The input FASTQ file and the
 * read ID strings are expected to have the same sorting. Entries in the read ID
 * strings may appear multiple times on consecutive lines.
 *
 * If the input filename ends in .gz, the input and output will be treated as
 * gzip-compressed (but not the read-ID list, which is read from STDIN).
 *
 *
 * This program was initially a PERL script, but this C++ version may give
 * slightly better performance.
 *
 *
 */

#define MAX_LEN 1024
#define MAX_QUEUE_SIZE 100

using namespace std;

mutex mx;
queue<vector<string>> q;
std::condition_variable cv;
bool finished = false, output_error = false;

void writerThread(ostream* output) {
    bool keep_running = true;
    while (keep_running) {
        vector<vector<string>> liness;
        {
            std::unique_lock<std::mutex> lk(mx);
            if (!finished)
                cv.wait(lk, []{return finished || !q.empty();});
            if (!q.empty()) {
                keep_running = true;
                do {
                    liness.emplace_back(move(q.front()));
                    q.pop();
                }
                while (!q.empty());
            }
            else {
                keep_running = !finished;
            }
            cv.notify_all();
        }
        for_each(liness.begin(), liness.end(), [&](vector<string> lines) {
            for_each(lines.begin(), lines.end(), [&](string s) { (*output) << s << '\n'; });
            }); // sorry
        if (!output) {
            output_error = true;
            return;
        }
    }
}


int main(int argc, char* argv[]) {

    // Main function: Reads arguments and calls analysisLoop
    if (argc != 2) {
        cerr << "usage: " << argv[0] << " INPUT_FASTQ_FILE" << endl;
        return 1;
    }
    string filename(argv[1]);

    istream* input_ptr;
    ostream* output_ptr;
    boost::iostreams::filtering_istream in;
    boost::iostreams::filtering_ostream out;

    bool valid;

    // Disable sync with printf, etc.
    ios_base::sync_with_stdio(false);
    // Use a simple locale, for speed
    setlocale(LC_ALL,"C");
    // Prevents flushing cout when reading from cin
    cin.tie(nullptr);

    ifstream file_input;
    file_input.open(filename, ios_base::in | ios_base::binary);
    if (file_input.fail()) {
        cerr << "Input file read error" << endl;
        return 1;
    }
    // Enable compression for input and output if file ends in ".gz"
    const string ending(".gz");
    if (equal(ending.rbegin(), ending.rend(), filename.rbegin())) {
        in.push(boost::iostreams::gzip_decompressor());
        in.push(file_input);
        input_ptr = &in;
        // Herw we set a gzip level of 4, used on many fastq files, and a 1 MB buffer
        out.push(boost::iostreams::gzip_compressor(4, 1024*1024));
        out.push(cout);
        output_ptr = &out;
    }
    else {
        input_ptr = &file_input;
        output_ptr = &cout;
    }
    istream & input = *input_ptr;

    thread outputer(writerThread, output_ptr);

    string data, header_tag;
    bool input_eof = false;
    while (!output_error && !input_eof) {
        bool accept_all;
        string old_header_tag(header_tag);
        // Remove repeated ID strings (these do happen in suprDUPr.read_id)
        while (header_tag == old_header_tag && !accept_all) {
            if (getline(cin, data)) {
                accept_all = false;
                size_t tab = data.find('\t');
                header_tag = data.substr(0, tab);
            }
            else if (cin.eof()) {
                accept_all = true;
            }
            else {
                cerr << "Input error while reading IDs from standard input" << endl;
                return 1;
            }
        }
        bool skippable_read_found = false;
        while (!output_error && !skippable_read_found && !input_eof) {
            vector<string> lines(4);
            for (int i=0; i<4; ++i) {
                getline(input, lines[i]);
                if (input.eof() && (lines[0] == "" && i==0)) {
                    input_eof = true;
                    break;
                }
                else if (!input) {
                    cerr << "error: Unexpected end of file while reading FASTQ file." << endl;
                    cerr << strerror(errno) << endl;
                    return 1;
                }
            }
            if (!input_eof) {
                // Check if header matches one of the skippable IDs from stdin
                if (!accept_all) {
                    string test = lines[0].substr(1, lines[0].find(' ')-1);
                    if (header_tag == test) {
                        skippable_read_found = true;
                    }
                }
                if (!skippable_read_found) {
                    unique_lock<mutex> lk(mx);
                    if (q.size() > MAX_QUEUE_SIZE) {
                        cv.wait(lk, [] {return q.size() < MAX_QUEUE_SIZE;});
                    }
                    q.push(move(lines));
                    cv.notify_all();
                }
            }
        }
    }
    {
        lock_guard<mutex> lk(mx);
        finished = true;
        cv.notify_all();
    }
    outputer.join();

    return 0;
}

