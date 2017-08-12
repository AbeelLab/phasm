#ifndef PHASM_OVERLAPPER_H
#define PHASM_OVERLAPPER_H

#include <string>
#include <iterator>
#include <vector>
#include <tuple>

#include <seqan/sequence.h>

using seqan::CharString;
using seqan::StringSet;
using std::string;
using std::vector;
using std::tuple;

typedef tuple<string, string, int, int, int, int> OverlapT;

class ExactOverlapper {
    public:
        ExactOverlapper();

        void addSequence(const string& id, const string& seq);
        vector<OverlapT> overlaps(unsigned int min_length);
    private:
        StringSet<CharString> readset;
        vector<string> readIds;
        vector<unsigned int> readLengths;
};

#endif
