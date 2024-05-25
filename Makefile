PROJECT_SOURCE_DIR ?= $(abspath ./)
PROJECT_NAME ?= $(shell basename $(PROJECT_SOURCE_DIR))
NUM_JOBS ?= 8

all:
	@echo nothing special

clean:
	rm -rf build dist src site wheelhouse *.egg-info
force_clean:
	docker run --rm -v `pwd`:`pwd` -w `pwd` -it alpine/make make clean
.PHONY: clean force_clean

package:
	mkdir -p cubao_headers
	cp -rf include cubao_headers
	python3 build.py > setup.py
