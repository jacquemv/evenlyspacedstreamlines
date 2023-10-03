OPTIONS = -Ofast -o bench -fopenmp -I../../evenlyspacedstreamlines -I../unittests

streamlines:
	g++ bench_streamlines.cpp ${OPTIONS} -DTIMING
	./bench > out

profile:
	g++ bench_streamlines.cpp ${OPTIONS} -pg -fno-pie -fno-builtin 
	./bench
	gprof ./bench > profile.txt
	rm bench gmon.out

triangularmesh:
	g++ bench_triangularmesh.cpp ${OPTIONS}
	./bench

maxsubinterval:
	g++ bench_maxsubinterval.cpp ${OPTIONS}
	./bench

memoryleaks:
	PYTHONPATH=../../:${PYTHONPATH} python test_memoryleaks.py

clean:
	rm -f out profile.txt bench
