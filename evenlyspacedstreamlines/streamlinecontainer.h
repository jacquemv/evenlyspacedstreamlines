#ifndef STREAMLINECONTAINER_H_
#define STREAMLINECONTAINER_H_

class StreamlineContainer {
public:
    int n;
    double* points;
    int* idx;

    StreamlineContainer();
    StreamlineContainer(int size);
    ~StreamlineContainer();

    void push_segment(double* P, double* Q, int face);
    int first_vector(double* vector);
    void reverse();
    double length();
    
    void print();

private:
    int nmax;
};

#endif