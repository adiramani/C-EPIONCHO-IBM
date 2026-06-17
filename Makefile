CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native
DEBUG_FLAGS = -std=c++17 -g -Wall -Wextra

SRC = $(wildcard epioncho_ibm/*.cpp)
BIN = c-epioncho-ibm

ifeq ($(OPENMP),1)

    CXX = /opt/homebrew/opt/llvm/bin/clang++

    OMP_INC = /opt/homebrew/opt/libomp/include
    OMP_LIB = /opt/homebrew/opt/libomp/lib

    CXXFLAGS += -fopenmp -isysroot $(shell xcrun --show-sdk-path)
    DEBUG_FLAGS += -fopenmp

    CPPFLAGS += -I$(OMP_INC)
    LDFLAGS += -L$(OMP_LIB) -lomp

endif

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

debug:
	$(CXX) $(DEBUG_FLAGS) $(SRC) -o $(BIN)-debug

clean:
	rm -f $(BIN) $(BIN)-debug

.PHONY: all debug clean
