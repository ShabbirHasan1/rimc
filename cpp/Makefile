CXX ?= clang++
CXXFLAGS = -std=c++11 -O3 -Wall
TARGET = cppsimc 
SRC = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

.PHONY: all clean
