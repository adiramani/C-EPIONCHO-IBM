SRC = $(wildcard epioncho_ibm/*.cpp)
BIN_DIR ?= .
BIN_NAME ?= c-epioncho-ibm
BIN = $(BIN_DIR)/$(BIN_NAME)

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
    
    CPU_MODEL := $(shell grep "model name" /proc/cpuinfo | head -1 | grep -o "EPYC\|Xeon")
    ifeq ($(CPU_MODEL),Xeon)
        # Ice Lake Xeon
        MARCH_FLAG = -march=icelake-server
    else ifeq ($(CPU_MODEL),EPYC)
        # AMD EPYC (assume Zen3 for 7742)
        MARCH_FLAG = -march=znver3
    else
        # Fallback
        MARCH_FLAG = -march=native
    endif

    CXXFLAGS = -std=c++17 -O3 $(MARCH_FLAG) -fopenmp -flto=auto
    DEBUG_FLAGS = -std=c++17 -g -Wall -Wextra -fopenmp

	GCC_MODULE := GCC/15.2.0

    SHELL := /bin/bash
    .SHELLFLAGS := -ec
    $(shell module load $(GCC_MODULE))
endif

all: $(BIN)

$(BIN): $(SRC)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

debug:
	mkdir -p $(BIN_DIR)
	$(CXX) $(DEBUG_FLAGS) $(CXXFLAGS) $(SRC) -o $(BIN)-debug

clean:
	rm -f $(BIN) $(BIN)-debug

.PHONY: all debug clean
