CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++17 -Iinclude

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
BIN_DIR = bin

TARGET = $(BIN_DIR)/main2

SRCS = $(SRC_DIR)/main2.cpp

OBJS = $(BUILD_DIR)/main2.o

all: $(TARGET)

$(TARGET): $(OBJS)
	@mkdir $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

.PHONY: all clean
