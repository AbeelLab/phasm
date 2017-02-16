#include <string>
#include <tuple>
#include <seqan/sequence.h>

#ifndef PHASM_ALIGNMENTS_H
#define PHASM_ALIGNMENTS_H

namespace phasm {

typedef std::vector<std::tuple<int, int> > TracePoints;

enum class Strand : uint8_t {
    SAME = 0,
    OPPOSITE = 1
};

/**
 * Represents a single read.
 *
 * Does not necessarily contain the sequence, only the length is required.
 * @tparam T String type of a read, for example `seqan::DnaString`
 */
template<typename T>
class Read {
public:
    Read(const std::string& id, unsigned long length);
    Read(const std::string& id, const T& sequence);

    unsigned long getLength() const;

private:
    std::string id;
    T sequence;
    unsigned long length;
};

template<typename T>
Read<T>::Read(const std::string& id, unsigned long length) : id(id), length(length) {
}

template<typename T>
Read<T>::Read(const std::string& id, const T& sequence) : id(id), sequence(sequence) {
    this->length = seqan::length(this->sequence);
}

template<typename T>
unsigned long Read<T>::getLength() const {
    return this->length;
}



/**
 * Represents a local alignment between two reads.
 *
 * Optionally includes trace points which can be used to reconstruct the actual
 * alignment.
 * @tparam T
 */
template<typename T>
class LocalAlignment {
public:
    LocalAlignment(const Read<T>& a, const Read<T>& b, Strand strand, const std::tuple<int, int>& arange,
                   const std::tuple<int, int>& brange);
    LocalAlignment(const Read<T>& a, const Read<T>& b, Strand strand, const std::tuple<int, int>& arange,
                   const std::tuple<int, int>& brange, const TracePoints& trace_points);

private:
    const Read<T>& a;
    const Read<T>& b;
    Strand strand;
    std::tuple<int, int> arange;
    std::tuple<int, int> brange;
    TracePoints trace_points;

};

template<typename T>
LocalAlignment<T>::LocalAlignment(const Read<T> &a, const Read<T> &b, Strand strand, const std::tuple<int, int>& arange,
                                  const std::tuple<int, int>& brange
)  : a(a), b(b), strand(strand), arange(arange), brange(brange)
{
}

template<typename T>
LocalAlignment<T>::LocalAlignment(const Read<T> &a, const Read<T> &b, Strand strand, const std::tuple<int, int>& arange,
                                  const std::tuple<int, int>& brange, const TracePoints& trace_points
)  : a(a), b(b), strand(strand), arange(arange), brange(brange), trace_points(trace_points)
{
}

} // namespace phasm

#endif //PHASM_ALIGNMENTS_H_H
