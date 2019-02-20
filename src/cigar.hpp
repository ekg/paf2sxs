#ifndef CIGAR_HPP_INCLUDED
#define CIGAR_HPP_INCLUDED

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

namespace paf2sxs {

struct cigar_op_t { uint64_t len; char op; };
typedef std::vector<cigar_op_t> cigar_t;
cigar_t cigar_from_string(const std::string& s);
std::string cigar_to_string(const cigar_t& cigar);

}

#endif
