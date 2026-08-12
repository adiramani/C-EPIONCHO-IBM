SRC = $(wildcard epioncho_ibm/*.cpp)
BIN = c-epioncho-ibm

OSNAME := $(shell uname -s)

ifeq ($(OSNAME),Linux)
    PLATFORM ?= linux
else
    PLATFORM ?= macos
endif

ifeq ($(PLATFORM),macos)
    CXX = /opt/homebrew/opt/llvm/bin/clang++
    CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp -isysroot $(shell xcrun --show-sdk-path)
    DEBUG_FLAGS = -std=c++17 -g -Wall -Wextra -fopenmp -isysroot $(shell xcrun --show-sdk-path)
endif

ifeq ($(PLATFORM),linux)
    CXX = g++
    CXXFLAGS = -std=c++17 -O3 -march=skylake -fopenmp
    DEBUG_FLAGS = -std=c++17 -g -Wall -Wextra -fopenmp

	GCC_MODULE := GCC/15.2.0

    SHELL := /bin/bash
    .SHELLFLAGS := -ec
    $(shell module load $(GCC_MODULE))
endif

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

debug:
	$(CXX) $(DEBUG_FLAGS) $(CXXFLAGS) $(SRC) -o $(BIN)-debug

clean:
	rm -f $(BIN) $(BIN)-debug

.PHONY: all debug clean
