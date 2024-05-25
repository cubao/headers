PROJECT_SOURCE_DIR ?= $(abspath ./)
PROJECT_NAME ?= $(shell basename $(PROJECT_SOURCE_DIR))
NUM_JOBS ?= 8

all:
	@echo nothing special

clean:
	rm -rf build dist cubao_headers site wheelhouse *.egg-info
force_clean:
	docker run --rm -v `pwd`:`pwd` -w `pwd` -it alpine/make make clean
.PHONY: clean force_clean

package:
	mkdir -p cubao_headers
	cp -rf include cubao_headers
	python3 build.py > setup.py
	python3 setup.py sdist

upload:
	python3 -m pip install twine
	twine upload dist/cubao_headers-*.tar.gz -r pypi
