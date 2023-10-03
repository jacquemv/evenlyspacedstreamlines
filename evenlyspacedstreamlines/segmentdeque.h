#ifndef SEGMENTDEQUE_H_
#define SEGMENTDEQUE_H_

//-----------------------------------------------------------------------------
typedef struct QueueElement {
    int idx;
    double s;
} QueueElement;

//-----------------------------------------------------------------------------
class SegmentDeque {
public:
    SegmentDeque();
    SegmentDeque(int maxsize_); // generate a buffer with 'maxsize' elements
    ~SegmentDeque();
    void clear(); // empty the queue

    // return 0 if appending/prepending is possible, and 1 otherwise
    inline bool append_failed();
    inline bool prepend_failed();
     // returns 1 if successful, 0 otherwise
    inline bool append(int idx, double s);
    inline bool prepend(int idx, double s);
    inline bool peek_first(int& idx, double& s);
    inline bool peek_position(int pos, int& idx, double& s);
    inline bool rewind(); // start a loop, returns 0 if empty
     // read and move the cursor to the next element
     // returns 0 if the end is reached
    inline bool peek_next(int& idx, double& s);
    // returns true if the cursor is at the first or last element
    inline bool is_extremity(); 
    // additional data at both extremities allows storage of partial segments
    void set_data(double f0, double f1, double l0, double l1);
    inline bool get_data(double& d1, double& d2);
    inline bool get_data(int pos, double& d1, double& d2);
    int size(); // number of elements in the queue
    // select a subset of the elements
    // start is the index from the bottom of the stack
    // to select all, start=0, len=size()
    void select(int start, int len);
    // flip the stack
    void flip();

//private:
    int maxsize;
    int first, last, cursor;
    QueueElement *array;
    double data_first[2], data_last[2];
};

#endif