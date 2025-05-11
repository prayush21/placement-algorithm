# Compiler settings
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -g
INCLUDES = -I./include
LDFLAGS =

# Directories
SRC_DIR = src
BUILD_DIR = build
INCLUDE_DIR = include

# Source files and Object files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Header files
HDRS = $(wildcard $(INCLUDE_DIR)/*.h)

# Target executable
TARGET = placer_tool

# Default target
all: $(BUILD_DIR)/$(TARGET)

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Link executable
$(BUILD_DIR)/$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

# Clean build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Run with default settings
run: $(BUILD_DIR)/$(TARGET)
	./$(BUILD_DIR)/$(TARGET) -h

# Format code
format:
	clang-format -i $(SRC_DIR)/*.cpp $(INCLUDE_DIR)/*.h

# Test target (placeholder)
test:
	@echo "No tests available yet."

# Dependencies
$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp $(INCLUDE_DIR)/parser.h $(INCLUDE_DIR)/placer.h $(INCLUDE_DIR)/fm.h $(INCLUDE_DIR)/router.h $(INCLUDE_DIR)/evaluator.h $(INCLUDE_DIR)/visualizer.h $(INCLUDE_DIR)/global_router.h $(INCLUDE_DIR)/detailed_router.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/parser.o: $(SRC_DIR)/parser.cpp $(INCLUDE_DIR)/parser.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/structures.o: $(SRC_DIR)/structures.cpp $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/fm.o: $(SRC_DIR)/fm.cpp $(INCLUDE_DIR)/fm.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/placer.o: $(SRC_DIR)/placer.cpp $(INCLUDE_DIR)/placer.h $(INCLUDE_DIR)/fm.h $(INCLUDE_DIR)/evaluator.h $(INCLUDE_DIR)/visualizer.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/visualizer.o: $(SRC_DIR)/visualizer.cpp $(INCLUDE_DIR)/visualizer.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/router.o: $(SRC_DIR)/router.cpp $(INCLUDE_DIR)/router.h $(INCLUDE_DIR)/global_router.h $(INCLUDE_DIR)/detailed_router.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/evaluator.o: $(SRC_DIR)/evaluator.cpp $(INCLUDE_DIR)/evaluator.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/global_router.o: $(SRC_DIR)/global_router.cpp $(INCLUDE_DIR)/global_router.h $(INCLUDE_DIR)/router.h $(INCLUDE_DIR)/structures.h
$(BUILD_DIR)/detailed_router.o: $(SRC_DIR)/detailed_router.cpp $(INCLUDE_DIR)/detailed_router.h $(INCLUDE_DIR)/router.h $(INCLUDE_DIR)/structures.h

.PHONY: all clean run format test 