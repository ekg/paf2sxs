#include <iostream>
#include <cassert>
#include "gzstream.h"
#include "args.hxx"
#include "paf.hpp"

int main(int argc, char** argv) {

    args::ArgumentParser parser("convert between PAF and SXS format, reading (possibly gzipped) PAF and writing .sxs on stdout");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> paf_filename(parser, "FILE", "read this PAF file (- for stdin)", {'p', "paf-in"});
    //args::Flag kmers_stdout(parser, "", "write the kmers to stdout", {'c', "stdout"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    std::string paf_file = args::get(paf_filename);
    if (paf_file == "-") paf_file = "/dev/stdin";
    igzstream paf_in(paf_file.c_str());
    if (!paf_in.good()) assert("PAF is not good!");
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(paf_in, line)) {
        paf2sxs::paf_row_t paf(line);
        // switch to sxs start/end format
        if (!paf.query_target_same_strand) {
            // we are reporting coordinates in the reverse complement
            // put them on the forward strand of the sequence
            // then flip them
            paf.query_start = paf.query_sequence_length - paf.query_start;
            paf.query_end = paf.query_sequence_length - paf.query_end;
            std::swap(paf.query_start, paf.query_end);
        }
        std::cout << "A" << "\t" << paf.target_sequence_name << "\t" << paf.query_sequence_name << "\n"
                  << "I" << "\t" << paf.target_start << "\t" << paf.target_end << "\t" << paf.query_start << "\t" << paf.query_end << "\n"
                  << "M" << "\t" << paf.num_matches << "\n"
                  << "C" << "\t" << paf2sxs::cigar_to_string(paf.cigar) << "\n"
                  << "Q" << "\t" << paf.mapping_quality << "\n";
    }
    return 0;
}
