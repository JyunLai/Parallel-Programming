appname := OMicroN

# 1. 編譯器設定
CXX := g++

# 2. 通用參數 (Common Flags)
# 注意：移除了 -fopenmp，我們稍後會針對不同版本分別加入
# 移除了 -mno-avx 以適應不同工作站架構
COMMON_FLAGS := -std=c++17 -ffast-math -Wno-unknown-pragmas -fPIC -fmax-include-depth=300

# 3. 路徑設定 (維持您之前的設定)
# HDF5, LibTIFF, Eigen 的路徑
PATH_HDF5 := $(HOME)/OMicroN/hdf5_local
PATH_TIFF := $(HOME)/OMicroN/libtiff_local
PATH_EIGEN := .

# Include 路徑
INCLUDES := -Isrc -Isrc/utils -I$(PATH_EIGEN) -I$(PATH_HDF5)/include -I$(PATH_TIFF)/include

# Library 路徑與連結
# 注意：Linux 下通常需要 -lhdf5_cpp 在 -lhdf5 之前
LIBS_DIRS := -L$(PATH_HDF5)/lib -L$(PATH_TIFF)/lib
LIBS_NAMES := -lm -ltiff -lpthread -lhdf5_cpp -lhdf5 -ldl

# 4. 原始碼設定
SRC_DIR := src
UTILS_DIR := $(SRC_DIR)/utils
SOURCES := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(UTILS_DIR)/*.cpp)
SOURCES_TMP := $(patsubst $(SRC_DIR)/%,%,$(SOURCES))
SOURCES_TMP := $(patsubst $(UTILS_DIR)/%,%,$(SOURCES_TMP))

# ==========================================
#        設定兩個版本的建置目標
# ==========================================

# --- Parallel Version (平行版) ---
BIN_PAR := bin/OMicroN_parallel
BUILD_PAR := build_parallel
FLAGS_PAR := $(COMMON_FLAGS) -fopenmp
LIBS_PAR  := $(LIBS_DIRS) $(LIBS_NAMES) -fopenmp
OBJS_PAR  := $(patsubst %.cpp,$(BUILD_PAR)/%.o,$(SOURCES_TMP))

# --- Serial Version (序列版/Baseline) ---
BIN_SER := bin/OMicroN_serial
BUILD_SER := build_serial
FLAGS_SER := $(COMMON_FLAGS) -Wno-unknown-pragmas # 忽略 #pragma omp
LIBS_SER  := $(LIBS_DIRS) $(LIBS_NAMES)           # 不連結 OpenMP
OBJS_SER  := $(patsubst %.cpp,$(BUILD_SER)/%.o,$(SOURCES_TMP))

# ==========================================
#             Make 指令規則
# ==========================================

# 預設執行 'make' 時，同時編譯兩個版本
all: parallel serial

parallel: $(BIN_PAR)
serial: $(BIN_SER)

# --- 建置 Parallel 執行檔 ---
$(BIN_PAR): $(OBJS_PAR) | bin
	$(CXX) -o $@ $^ $(LIBS_PAR)
	@echo "Build Success: Parallel version is at $@"

# --- 建置 Serial 執行檔 ---
$(BIN_SER): $(OBJS_SER) | bin
	$(CXX) -o $@ $^ $(LIBS_SER)
	@echo "Build Success: Serial version is at $@"

# --- 編譯 Parallel 物件檔 (.o) ---
$(BUILD_PAR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_PAR)
	@mkdir -p $(dir $@)
	$(CXX) -c $(FLAGS_PAR) $(INCLUDES) -o $@ $<

$(BUILD_PAR)/%.o: $(UTILS_DIR)/%.cpp | $(BUILD_PAR)
	@mkdir -p $(dir $@)
	$(CXX) -c $(FLAGS_PAR) $(INCLUDES) -o $@ $<

# --- 編譯 Serial 物件檔 (.o) ---
$(BUILD_SER)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_SER)
	@mkdir -p $(dir $@)
	$(CXX) -c $(FLAGS_SER) $(INCLUDES) -o $@ $<

$(BUILD_SER)/%.o: $(UTILS_DIR)/%.cpp | $(BUILD_SER)
	@mkdir -p $(dir $@)
	$(CXX) -c $(FLAGS_SER) $(INCLUDES) -o $@ $<

# --- 建立資料夾 ---
bin $(BUILD_PAR) $(BUILD_SER):
	mkdir -p $@

# --- 清除指令 ---
clean:
	rm -rf bin build_parallel build_serial

.PHONY: all parallel serial clean