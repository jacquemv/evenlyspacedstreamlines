#include "engine.cpp"

void test_basic()
{
    StreamlineContainer sc(4);
    for (int i=0;i<3;i++) {
        double P[3] = {0, 0, 0}, Q[3] = {0, 0, 0};
        P[0] = (i+1)*0.5;
        Q[0] = i*0.5;
        sc.push_segment(P, Q, i);
    }
    assert( sc.idx[0] == 0 && sc.idx[1] == 1 && sc.idx[2] == 2 );
    assert( sc.length() == 1.5 );
}

void test_connect()
{
    double a[] = {1, 2, 3};
    double b[] = {2, -2, 4};
    double c[] = {5, 2, -1};
    double d[] = {3, 1, 0};

    StreamlineContainer s1(4);
    s1.push_segment(a, b, 0);
    s1.push_segment(b, c, 0);
    s1.push_segment(c, d, 0);

    StreamlineContainer s2(4);
    s2.push_segment(b, a, 0);
    s2.push_segment(b, c, 0);
    s2.push_segment(c, d, 0);
    assert( s1.length() == s2.length() );

    StreamlineContainer s3(4);
    s3.push_segment(b, a, 0);
    s3.push_segment(c, b, 0);
    s3.push_segment(d, c, 0);
    assert( s1.length() == s3.length() );
}


int main()
{
    test_basic();
    test_connect();
}