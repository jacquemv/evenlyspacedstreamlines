OPTIONS = -I.. -O3 -o test -fopenmp -I../../evenlyspacedstreamlines

AUTOTESTS = connections disjointintervals quadraticforms segmentdeque \
			streamlinecontainer triangle maxsubinterval

auto: ${AUTOTESTS} triangularmesh segmentsearch wrapper triangle_idx
	rm -f ./test out

${AUTOTESTS}:
	g++ test_$@.cpp ${OPTIONS}
	./test

triangularmesh:
	g++ test_$@.cpp ${OPTIONS}
	./test
	python test_triangularmesh.py
	rm -f test_mesh3d.{adjtri,inctri,neightri,radius,bound,nbound}

segmentsearch:
	g++ test_$@.cpp ${OPTIONS}
	./test > out
	diff out test_segmentsearch_ref.txt

wrapper:
	PYTHONPATH=../../:${PYTHONPATH} python test_wrapper.py

triangle_idx:
	PYTHONPATH=../../:../visualtests:${PYTHONPATH} python test_triangle_idx.py