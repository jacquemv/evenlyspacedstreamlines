all:
	python -m build #--no-isolation

local:
	python setup.py build_ext -i
	mv *.so evenlyspacedstreamlines

clean:
	rm -rf build dist *.egg-info