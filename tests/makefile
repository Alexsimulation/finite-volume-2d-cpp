MPICC := mpic++

SOURCES := $(shell find $(FVHYPER_DIR)/fvhyper/src -name '*.cpp')
INCLUDES := -I${FVHYPER_DIR}
OPTIM := -O3

build:
	${MPICC} -o tests tests.cpp ${SOURCES} ${INCLUDES} ${OPTIM} -std=c++17
