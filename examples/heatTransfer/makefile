MPICC?=mpic++


FVHYPER_DIR := ../..
SOURCEDIR := ${FVHYPER_DIR}/fvhyper/src
SOURCES := $(shell find $(SOURCEDIR) -name '*.cpp')

INCLUDES := -I${FVHYPER_DIR}
OPTIM := -O3

build:
	@echo ${SOURCES}
	${MPICC} -o main main.cpp ${SOURCES} ${INCLUDES} ${OPTIM} -std=c++17
