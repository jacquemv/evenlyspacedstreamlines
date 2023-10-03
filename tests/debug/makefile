OPTIONS = -O3 -o test -fopenmp -I../../evenlyspacedstreamlines -I../unittests

streamlines:
	g++ debug_streamlines.cpp ${OPTIONS}
	mkdir -p output
	./test
	rm test
	python plot_streamlines.py ../data/test_mesh2d.off 0.2

periodic:
	g++ debug_periodic.cpp ${OPTIONS}
	mkdir -p output
	./test
	rm test
	python plot_streamlines.py ../data/test_annulus.off 0.04

uturn:
	g++ debug_uturn.cpp ${OPTIONS}
	mkdir -p output
	./test
	rm test
	python plot_streamlines.py ../data/test_uturn.off 0.2

clean:
	rm -rf output
