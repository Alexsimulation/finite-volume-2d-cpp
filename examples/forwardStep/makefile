MPICC := mpic++

SOURCES := $(shell find $(FVHYPER_DIR)/fvhyper/src -name '*.cpp')
INCLUDES := -I${FVHYPER_DIR}
OPTIM := -O3 -march=native -fno-math-errno

build:
	${MPICC} -o main main.cpp ${SOURCES} ${INCLUDES} ${OPTIM} -std=c++17
