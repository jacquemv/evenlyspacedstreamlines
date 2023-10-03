OPTIONS = -I.. -O3 -o test -fopenmp -I../../evenlyspacedstreamlines -I../unittests

all: curvature segmentsearch streamlines periodic klein neighborhood

neighborhood:
	g++ test_neighborhood.cpp ${OPTIONS}
	./test
	python test_neighborhood.py

curvature:
	g++ test_curvature.cpp ${OPTIONS}
	./test > out
	python test_curvature.py
	rm out test

segmentsearch:
	g++ test_segmentsearch.cpp ${OPTIONS}
	python test_segmentsearch.py
	rm test

streamlines:
	g++ test_streamlines.cpp ${OPTIONS}
	./test > out
	python test_streamlines.py
	rm test out

periodic:
	PYTHONPATH=../../:${PYTHONPATH} python test_periodic.py

klein:
	PYTHONPATH=../../:${PYTHONPATH} python test_klein.py