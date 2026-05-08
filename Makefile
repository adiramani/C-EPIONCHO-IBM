CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ljsoncpp
DEBUG_FLAGS = -std=c++17 -g -Wall -Wextra

SRC = $(wildcard epioncho_ibm/*.cpp)
BIN = c-epioncho-ibm

all: $(BIN)

$(BIN): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

debug:
	$(CXX) $(DEBUG_FLAGS) $(SRC) -o $(BIN)-debug

clean:
	rm -f $(BIN) $(BIN)-debug

.PHONY: all debug clean
