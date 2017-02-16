#include <string>
#include <tuple>
#include <algorithm>
#include <seqan/sequence.h>

#ifndef PHASM_ALIGNMENTS_H
#define PHASM_ALIGNMENTS_H

namespace phasm
{

using std::get;
using seqan::length;

typedef std::vector<std::tuple<int, int> > TracePoints;

/**
 * Distinguish if two reads come from the same or opposite strand.
 */
enum class Strand : uint8_t
{
    SAME,
    OPPOSITE
};


/**
 * Distinguish between different kind of overlaps.
 *
 * An overlap is contained if one read completely spans the other.
 */
enum class AlignmentType : uint8_t
{
    A_CONTAINED,
    B_CONTAINED,
    OVERLAP_AB,
    OVERLAP_BA
};

/**
 * Represents a single read.
 *
 * Does not necessarily contain the sequence, only the length is required.
 * @tparam T String type of a read, for example `seqan::DnaString`
 */
template<typename T>
class Read
{
public:
    Read(const std::string &id, unsigned long length);
    Read(const std::string &id, const T &sequence);

    unsigned long length() const;

private:
    std::string id;
    T sequence;
    unsigned long length_;
};

template<typename T>
Read<T>::Read(const std::string &id, unsigned long length) : id(id), length_(length)
{
}

template<typename T>
Read<T>::Read(const std::string &id, const T &sequence) : id(id), sequence(sequence)
{
    length_ = seqan::length(sequence);
}

template<typename T>
unsigned long Read<T>::length() const
{
    return length_;
}


/**
 * Represents a local alignment between two reads.
 *
 * Optionally includes trace points which can be used to reconstruct the actual
 * alignment.
 * @tparam T
 */
template<typename T>
class LocalAlignment
{
public:
    LocalAlignment(const Read<T> &a, const Read<T> &b, Strand strand, const std::tuple<int, int> &arange,
                   const std::tuple<int, int> &brange, int differences);

    LocalAlignment(const Read<T> &a, const Read<T> &b, Strand strand, const std::tuple<int, int> &arange,
                   const std::tuple<int, int> &brange, int differences, const TracePoints &trace_points);

    const Read<T> &readA() const;
    const Read<T> &readB() const;

    const std::tuple<int, int> &getArange() const;
    const std::tuple<int, int> &getBrange() const;

    int getOverlapLength() const;
    int getDifferences() const;
    float getErrorRate() const;
    unsigned long getOverhang() const;

    AlignmentType classify() const;

private:
    const Read<T> &a;
    const Read<T> &b;
    Strand strand;
    std::tuple<int, int> arange;
    std::tuple<int, int> brange;
    int differences;
    TracePoints trace_points;

};

template<typename T>
LocalAlignment<T>::LocalAlignment(const Read<T> &a, const Read<T> &b, Strand strand, const std::tuple<int, int> &arange,
                                  const std::tuple<int, int> &brange, int differences
)  : a(a), b(b), strand(strand), arange(arange), brange(brange), differences(differences)
{
}

template<typename T>
LocalAlignment<T>::LocalAlignment(const Read<T> &a, const Read<T> &b, Strand strand, const std::tuple<int, int> &arange,
                                  const std::tuple<int, int> &brange, int differences, const TracePoints &trace_points
)  : a(a), b(b), strand(strand), arange(arange), brange(brange), differences(differences), trace_points(trace_points)
{
}

template<typename T>
const Read<T> &LocalAlignment<T>::readA() const
{
    return a;
}

template<typename T>
const Read<T> &LocalAlignment<T>::readB() const
{
    return b;
}

template<typename T>
const std::tuple<int, int> &LocalAlignment<T>::getArange() const
{
    return arange;
}

template<typename T>
const std::tuple<int, int> &LocalAlignment<T>::getBrange() const
{
    return brange;
}

template<typename T>
int LocalAlignment<T>::getOverlapLength() const
{
    return std::min(get<1>(arange) - get<0>(arange), get<1>(brange) - get<0>(brange));
}

template<typename T>
unsigned long LocalAlignment<T>::getOverhang() const
{
    return (std::min(get<0>(arange), get<0>(brange)) +
            std::min(a.length() - get<1>(arange), b.length() - get<1>(brange)));
}

template<typename T>
AlignmentType LocalAlignment<T>::classify() const
{
    if (get<0>(arange) <= get<0>(brange) && (a.length() - get<1>(arange)) <= (b.length() - get<1>(brange)))
    {
        return AlignmentType::A_CONTAINED;
    }
    else if (get<0>(arange) >= get<0>(brange) &&
             (a.length() - get<1>(arange)) >= (b.length() - get<1>(brange)))
    {
        return AlignmentType::B_CONTAINED;
    }
    else if (get<0>(arange) > get<0>(brange))
    {
        return AlignmentType::OVERLAP_AB;
    }
    else
    {
        return AlignmentType::OVERLAP_BA;
    }
}


} // namespace phasm

#endif //PHASM_ALIGNMENTS_H_H
