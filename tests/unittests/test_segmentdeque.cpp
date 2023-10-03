#include "engine.cpp"

double sumprod(SegmentDeque* queue)
{
    int idx;
    double s, x = 0;
    queue->rewind();
    while (queue->peek_next(idx, s))
        x += idx * s;
    return x;
}

#undef QUEUE_FOR_EACH_ELEMENT
#define QUEUE_FOR_EACH_ELEMENT(queue, idx, s) \
    queue.rewind(); \
    int idx; double s; \
    while (queue.peek_next(idx, s))

int main()
{
    int idx;
    double s;
    SegmentDeque queue(5);
    assert( queue.peek_first(idx, s) == 0);
    assert( queue.rewind() == 0);
    assert( queue.size() == 0 );

    assert( queue.append_failed() == 0 );
    assert( queue.append(0, 0.0) == 1);
    assert( queue.size() == 1 );
    assert( queue.append(1, 0.1) == 1);
    assert( queue.append_failed() == 0 );
    assert( queue.size() == 2 );
    assert( queue.append(2, 0.2) == 1);
    assert( queue.append_failed() == 1 );
    assert( queue.append(3, 0.3) == 0);
    assert( queue.prepend(-1, -0.1) == 1);
    assert( queue.prepend_failed() == 0 );
    assert( queue.prepend(-2, -0.2) == 1);
    assert( queue.prepend_failed() == 1 );
    assert( queue.prepend(-3, -0.3) == 0);
    assert( queue.prepend_failed() == 1 );
    assert( queue.size() == 5 );

    assert( queue.rewind() == 1 );
    int i = -2;
    while (queue.peek_next(idx, s)) {
        assert( (i == idx) || (i/10 == s) );
        assert( queue.is_extremity() == ( i == -2 || i == 2) );
        i++;
    }

    queue.peek_position(0, idx, s);  assert( idx == -2 && s == -0.2 );
    queue.peek_position(1, idx, s);  assert( idx == -1 && s == -0.1 );
    queue.peek_position(-1, idx, s);  assert( idx == 2 && s == 0.2 );
    queue.peek_position(-2, idx, s);  assert( idx == 1 && s == 0.1 );
    assert( queue.peek_position(5, idx, s) == 0 );

    double x = sumprod(&queue);
    queue.data_first[0] = 1;
    queue.data_first[1] = 2;
    queue.data_last[0] = 3;
    queue.data_last[1] = 4;
    int size_ref = queue.size();
    queue.flip();
    assert( queue.size() == size_ref );
    assert( sumprod(&queue) == x );
    assert( queue.data_first[0] == 3 );
    assert( queue.data_first[1] == 4 );
    assert( queue.data_last[0] == 1 );
    assert( queue.data_last[1] == 2 );

    x = 0;
    QUEUE_FOR_EACH_ELEMENT(queue, idx2, s2) {
        x += idx2 * s2;
        double d1, d2;
        if (queue.is_extremity()) {
            assert( queue.get_data(d1, d2) == 1 );
            if (idx2 == 2)
                assert( d1 == 3 && d2 == 4 );
            else if (idx2 == -2)
                assert( d1 == 1 && d2 == 2 );
            else assert(0);
        } else {
            assert( queue.get_data(d1, d2) == 0 );
        }
    }
    assert( sumprod(&queue) == x );

    queue.clear();
    assert( queue.rewind() == 0);
    assert( queue.peek_first(idx, s) == 0);

    assert( queue.append(0, 0.0) == 1);
    assert( queue.append(1, 0.1) == 1);
    assert( queue.prepend(-1, -0.1) == 1);
    assert( queue.size() == 3 );

    i = 0;
    assert( queue.rewind() == 1 );
    while (queue.peek_next(idx, s)) i++;
    assert( i == 3 );

    SegmentDeque q(128);
    for (int k=0;k<20;k++)
        q.append(k, 0.0);
    assert( q.size() == 20 );
    q.select(4, 8);
    assert( q.size() == 8 );
    assert( q.peek_first(idx, s) == 1);
    assert( idx == 4 );
    q.select(3, 5);
    assert( q.size() == 5 );
    assert( q.peek_first(idx, s) == 1);
    assert( idx == 7 );
    q.select(0, q.size());
    assert( q.peek_first(idx, s) == 1);
    assert( idx == 7 );
}