# 定义编译器
CXX = g++
#CXXFLAGS = -Wall -g -std=c++17  -DDEBUG# 编译选项，可以根据需要调整
CXXFLAGS = -Wall -O3 -std=c++17
# 定义目录
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

# 定义源文件和目标文件
SRCS = $(SRC_DIR)/finta.cpp $(SRC_DIR)/finta_fast.cpp $(SRC_DIR)/quantum_chip.cpp $(SRC_DIR)/circuit.cpp $(SRC_DIR)/token_swap.cpp
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# 定义最终的可执行文件名
FINTA = $(BIN_DIR)/finta
FINTA_FAST = $(BIN_DIR)/finta_fast

# 默认目标
all: $(FINTA) $(FINTA_FAST)

# 编译规则：将每个 .cpp 文件编译成对应的 .o 文件
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 创建中间文件目录
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# 链接规则：将所有的 .o 文件链接成最终的可执行文件
$(FINTA): $(BUILD_DIR)/finta.o $(BUILD_DIR)/quantum_chip.o $(BUILD_DIR)/circuit.o | $(BIN_DIR)
	$(CXX) $(BUILD_DIR)/finta.o $(BUILD_DIR)/quantum_chip.o $(BUILD_DIR)/circuit.o -o $@

$(FINTA_FAST): $(BUILD_DIR)/finta_fast.o $(BUILD_DIR)/quantum_chip.o $(BUILD_DIR)/circuit.o | $(BIN_DIR)
	$(CXX) $(BUILD_DIR)/finta_fast.o $(BUILD_DIR)/quantum_chip.o $(BUILD_DIR)/circuit.o -o $@

# 创建可执行文件目录
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# 清理规则：清除所有生成的 .o 文件和可执行文件
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# 目标依赖
$(BUILD_DIR)/quantum_chip.o: $(SRC_DIR)/quantum_chip.cpp
$(BUILD_DIR)/circuit.o: $(SRC_DIR)/circuit.cpp
$(BUILD_DIR)/finta.o: $(SRC_DIR)/finta.cpp
$(BUILD_DIR)/finta_fast.o: $(SRC_DIR)/finta_fast.cpp
$(BUILD_DIR)/token_swap.o: $(SRC_DIR)/token_swap.cpp
.PHONY: all clean
