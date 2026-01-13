.PHONY: all miniQuant clean

CXX := g++
CXXFLAGS := -g -fno-lto -std=c++23 -O3
LDFLAGS := -lpthread -ldl -lbz2
INCLUDE_DIR := ./include
LIB_DIR := ./lib
LIBS := -lggcat_api -lggcat_cpp_bindings -lggcat_cxx_interop -lz

SRC_DIR := ./src
SOURCES := main.cc util.cc align.cc calKvalue.cc EM.cc sequenceReader.cc alignmentReader.cc annotationParser.cc regionGenerator.cc condProbEstimator.cc lengthDensityEstimator.cc fileWriter.cc
OBJECTS := $(addprefix $(SRC_DIR)/, $(SOURCES:.cc=.o))

OUTPUT := miniQuant

all: miniQuant

miniQuant: $(OBJECTS) $(LIB_DIR)/libggcat_cpp_bindings.a $(LIB_DIR)/libggcat_cxx_interop.a $(LIB_DIR)/libggcat_api.a $(LIB_DIR)/libz.a
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -L$(LIB_DIR) $(OBJECTS) $(LIBS) -o $(OUTPUT) $(LDFLAGS)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $< -o $@

clean:
	rm -f $(OUTPUT) $(OBJECTS)