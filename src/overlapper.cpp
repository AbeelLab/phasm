#include <algorithm>
#include <unordered_map>
#include <stack>
#include <seqan/index.h>

#include "overlapper.h"

using std::unordered_map;
using std::stack;

using seqan::String;
using seqan::Index;
using seqan::IndexEsa;
using seqan::Iterator;
using seqan::TopDown;
using seqan::ParentLinks;
using seqan::EmptyEdges;

ExactOverlapper::ExactOverlapper() { };


void ExactOverlapper::addSequence(const string& id, const string& seq) {
    appendValue(this->readset, seq);
    this->readIds.push_back(id);
    this->readLengths.push_back(seq.length());
}

vector<OverlapT> ExactOverlapper::overlaps(unsigned int min_length) {
    vector<OverlapT> overlaps;
    unordered_map<int, stack<unsigned int> > read_stacks;

    // Create index using enhanced suffix array
    typedef Index<StringSet<CharString>, IndexEsa<> > TIndex;
    TIndex read_index(this->readset);
    typedef Iterator<TIndex, TopDown<ParentLinks<EmptyEdges> > >::Type TIterator;
    TIterator it(read_index);

    do {
        auto occ = getOccurrence(it);
        if(!isLeaf(it) && repLength(it) >= min_length) {
            // Inspect children, and check for empty edges (denotes an actual
            // suffix)
            goDown(it);
            do {
                if(isLeaf(it) && emptyParentEdge(it)) {
                    // We found an actual suffix of one of the strings,
                    // determine which one
                    auto occ = getOccurrence(it);
                    auto stack_it = read_stacks.find(occ.i1);
                    if(stack_it == read_stacks.end()) {
                        // No stack exists yet
                        read_stacks[occ.i1] = stack<unsigned int>();
                        read_stacks[occ.i1].push(repLength(it));
                    } else {
                        stack_it->second.push(repLength(it));
                    }
                }
            } while(goRight(it));

            // All children checked, go back to parent
            goUp(it);
        }

        if(isLeaf(it) && occ.i2 == 0) {
            // This is a leaf corresponding to one of the complete original
            // strings.
            // Each non-empty stack is a suffix-prefix overlap
            std::for_each(read_stacks.begin(), read_stacks.end(),
                    [&occ, &overlaps, this](auto& e) {
                        unsigned int other_id = e.first;

                        if(other_id == occ.i1) {
                            return;
                        }

                        stack<unsigned int>& read_stack = e.second;
                        auto suffix_len = read_stack.top();

                        unsigned int astart = this->readLengths[other_id] - suffix_len;
                        unsigned int aend = this->readLengths[other_id];
                        unsigned int bstart = 0;
                        unsigned int bend = suffix_len;

                        auto overlap = std::make_tuple(
                                this->readIds[other_id], this->readIds[occ.i1],
                                astart, aend, bstart, bend
                        );

                        overlaps.push_back(overlap);
                    }
            );

            // Check if this leaf is contained in another string, and also
            // output contained overlaps
            if(emptyParentEdge(it)) {
                // Visit siblings
                auto parent_it = TIterator(it);
                goUp(parent_it);
                auto contained_in = getOccurrences(parent_it);
                for(unsigned int i = 0; i < length(contained_in); ++i) {
                    auto other_id = contained_in[i].i1;
                    if(other_id == occ.i1) {
                        continue;
                    }

                    auto start_pos = contained_in[i].i2;
                    auto overlap = std::make_tuple(
                            this->readIds[other_id], this->readIds[occ.i1],
                            start_pos, start_pos + repLength(parent_it),
                            0, repLength(parent_it)
                    );

                    overlaps.push_back(overlap);
                }
            }
        }

        // Actual tree traversal part
        // Pre-order, but when visiting a node for the last time,
        // clean up read stacks.
        if (!goDown(it) && !goRight(it)) {
            do {
                // Last time we will visit the parent node,
                // check if we need to pop from read stacks.
                //
                // There's no goLeft, so goUp and goDown again
                goUp(it);
                goDown(it);
                do {
                    if(isLeaf(it) && emptyParentEdge(it)) {
                        // We found an actual suffix of one of the strings,
                        // determine which one
                        auto occ = getOccurrence(it);
                        if(read_stacks.find(occ.i1) != read_stacks.end()) {
                            read_stacks[occ.i1].pop();
                            if(read_stacks[occ.i1].empty()) {
                                read_stacks.erase(occ.i1);
                            }
                        }
                    }
                } while(goRight(it));

                // Go back up again
                goUp(it);
            } while (!goRight(it) && !isRoot(it));
        }
    } while(!isRoot(it));

    return overlaps;
}

