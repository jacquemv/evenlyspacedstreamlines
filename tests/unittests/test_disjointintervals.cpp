#include "engine.cpp"

const double EPS = 1e-10;

void test_reallocate()
{
    DisjointIntervals mask;
    for (int i=0;i<1000;i++) mask.insert(i/1000.0, (i+0.5)/1000.0);
    assert( mask.n == 1000 && mask.t[1998] == 0.999 );
}

void test_insert0(double a, double b, int n_expected)
{
    DisjointIntervals mask;
    mask.insert(a, b);
    assert( mask.n == n_expected );
}

void test_insert2(double a, double b, int n_expected)
{
    DisjointIntervals mask;
    mask.insert(0.1, 0.2);
    mask.insert(0.3, 0.4);
    mask.insert(a, b);
    assert( mask.n == n_expected );
}

void test_insert5(double a, double b, int n_expected)
{
    DisjointIntervals mask;
    mask.insert(0.1, 0.2);
    mask.insert(0.3, 0.4);
    mask.insert(0.5, 0.6);
    mask.insert(0.7, 0.8);
    mask.insert(0.9, 1.0);
    mask.insert(a, b);
    assert( mask.n == n_expected );
}

void test_insert_random()
{
    DisjointIntervals mask;
    for (int i=0;i<100;i++) {
        double a = ((double)rand()/(double)(RAND_MAX));
        double b = ((double)rand()/(double)(RAND_MAX));
        mask.insert(a, b);
        for (int j=1;j<2*mask.n;j++)
            assert( mask.t[j-1] < mask.t[j] );
    }
}

void test_contains0()
{
    DisjointIntervals mask;
    assert( !mask.contains(0.5) );
}

void test_contains1()
{
    DisjointIntervals mask;
    mask.insert(0.1, 0.2);
    for (int i=0;i<100;i++) {
        double x = i/100.0;
        assert( mask.contains(x) == ((x>=0.1 && x<=0.2)) );
    }
}

void test_contains2()
{
    DisjointIntervals mask;
    mask.insert(0.1, 0.2);
    mask.insert(0.3, 0.4);
    assert( mask.contains(0.2) );
    assert( mask.contains(0.3) );
    for (int i=0;i<99;i++) {
        double x = i/100.0+0.0005;
        assert( mask.contains(x) == ((x>=0.1 && x<=0.2) || (x>=0.3 && x<=0.4)) );
    }
}

void test_random()
{
    DisjointIntervals mask;
    for (int i=0;i<1000;i++) {
        double x = mask.pick_random();
        assert( x >= 0 && x <= 1 );
    }

    mask.insert(0.1, 0.2);
    mask.insert(0.3, 0.4);
    for (int i=0;i<1000;i++) {
        double x = mask.pick_random();
        assert( x >= 0 && x <= 1 );
        assert( !mask.contains(x) );
    }
}

void test_coverage()
{
    DisjointIntervals mask;
    assert( mask.is_empty() );
    assert( mask.coverage() == 0.0 );
    mask.insert(0.1, 0.2);
    assert( !mask.is_empty() );
    assert( fabs(mask.coverage() - 0.1) < EPS );
    mask.insert(0.3, 0.4);
    assert( fabs(mask.coverage() - 0.2) < EPS );
    mask.insert(0.35, 2);
    assert( fabs(mask.coverage() - 0.8) < EPS );
    mask.insert(-1, 0.5);
    assert( mask.coverage() == 1.0 );
    assert( mask.is_full() );
}


int main()
{
    test_reallocate();

    test_insert0(0.4, 0.7, 1);

    test_insert2(0.2, 0.25, 2);
    test_insert2(0.25, 0.3, 2);
    test_insert2(0.2, 0.3, 1);
    test_insert2(0.1, 0.4, 1);
    test_insert2(0.15, 0.31, 1);
    test_insert2(0.05, 0.6, 1);
    test_insert2(0.25, 0.6, 2);
    test_insert2(0.33, 0.34, 2);
    test_insert2(0.05, 0.06, 3);
    test_insert2(0.25, 0.26, 3);
    test_insert2(0.8, 0.9, 3);

    test_insert5(0.25, 0.85, 3);
    test_insert5(0.01, 0.75, 2);
    test_insert5(0.25, 2.8, 2);
    test_insert5(-5, 5, 1);

    test_contains0();
    test_contains1();
    test_contains2();
    
    test_insert_random();
    test_random();
    
    test_coverage();
}