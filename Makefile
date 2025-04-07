# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
INCLUDES = -I./include

# Directories
SRC_DIR = src
BUILD_DIR = build

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Target executable
TARGET = placer

# Default target
all: $(BUILD_DIR)/$(TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link executable
$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $@

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Run with default settings
run: $(BUILD_DIR)/$(TARGET)
	./$(BUILD_DIR)/$(TARGET) -h

# Dependencies
$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp include/parser.h include/placer.h include/fm.h
$(BUILD_DIR)/parser.o: $(SRC_DIR)/parser.cpp include/parser.h include/structures.h
$(BUILD_DIR)/structures.o: $(SRC_DIR)/structures.cpp include/structures.h
$(BUILD_DIR)/fm.o: $(SRC_DIR)/fm.cpp include/fm.h include/structures.h
$(BUILD_DIR)/placer.o: $(SRC_DIR)/placer.cpp include/placer.h include/fm.h include/structures.h
$(BUILD_DIR)/visualizer.o: $(SRC_DIR)/visualizer.cpp include/visualizer.h include/structures.h

.PHONY: all clean run 