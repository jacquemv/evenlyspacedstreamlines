#include "segmentdeque.h"

//-----------------------------------------------------------------------------
SegmentDeque::SegmentDeque()
{
    maxsize = 0;
    array = NULL;
}

//-----------------------------------------------------------------------------
SegmentDeque::SegmentDeque(int maxsize_)
{
    maxsize = maxsize_;
    array = new QueueElement [maxsize];
    clear();
}

//-----------------------------------------------------------------------------
SegmentDeque::~SegmentDeque()
{
    delete [] array;
}

//-----------------------------------------------------------------------------
void SegmentDeque::clear()
{
    first = cursor = maxsize / 2;
    last = first - 1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::append_failed()
{
    return last >= maxsize-1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::prepend_failed()
{
    return first <= 0;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::append(int idx, double s)
{
    if (last >= maxsize-1) return 0;
    last++;
    array[last].idx = idx;
    array[last].s = s;
    return 1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::prepend(int idx, double s)
{
    if (first <= 0) return 0;
    first--;
    array[first].idx = idx;
    array[first].s = s;
    return 1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::peek_first(int& idx, double& s)
{
    if (first > last) return 0;
    idx = array[first].idx;
    s = array[first].s;
    return 1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::peek_position(int pos, int& idx, double& s)
{
    if (first > last) return 0;
    if (pos < 0) pos += last + 1; else pos += first;
    if (pos < first || pos > last) return 0;
    idx = array[pos].idx;
    s = array[pos].s;
    return 1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::rewind()
{
    cursor = first;
    return last >= first;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::peek_next(int& idx, double& s)
{
    if (cursor > last) return 0;
    idx = array[cursor].idx;
    s = array[cursor].s;
    cursor++;
    return 1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::is_extremity()
{
    return cursor-1 == first || cursor-1 == last;
}

//-----------------------------------------------------------------------------
void SegmentDeque::set_data(double f0, double f1, double l0, double l1)
{
    data_first[0] = f0; data_first[1] = f1;
    data_last[0] = l0;  data_last[1] = l1;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::get_data(double& d1, double& d2)
{
    if (cursor-1 == first) {
        d1 = data_first[0];
        d2 = data_first[1];
        return 1;
    }
    if (cursor-1 == last) {
        d1 = data_last[0];
        d2 = data_last[1];
        return 1;
    }
    return 0;
}

//-----------------------------------------------------------------------------
inline bool SegmentDeque::get_data(int pos, double& d1, double& d2)
{
    if (pos == first) {
        d1 = data_first[0];
        d2 = data_first[1];
        return 1;
    }
    if (pos == last) {
        d1 = data_last[0];
        d2 = data_last[1];
        return 1;
    }
    return 0;
}

//-----------------------------------------------------------------------------
int SegmentDeque::size()
{
    return last - first + 1;
}

//-----------------------------------------------------------------------------
void SegmentDeque::select(int start, int len)
{
    first += start;
    last = first + len - 1;
    assert( first >= 0 && last < maxsize );
}

//-----------------------------------------------------------------------------
void SegmentDeque::flip()
{
    #define SWAP(x, y, tmp) tmp = x; x = y; y = tmp;
    // this function is only used when the array is at least half-full
    int i = 0, j = maxsize-1;
    while (i < j) {
        int tmpi;
        double tmpd;
        SWAP(array[i].idx, array[j].idx, tmpi);
        SWAP(array[i].s, array[j].s, tmpd);
        i++;
        j--;
    }
    int last_bkp = last;
    last = maxsize-1 - first;
    first = maxsize-1 - last_bkp;
    double tmpd;
    SWAP(data_first[0], data_last[0], tmpd);
    SWAP(data_first[1], data_last[1], tmpd);
    #undef SWAP
}