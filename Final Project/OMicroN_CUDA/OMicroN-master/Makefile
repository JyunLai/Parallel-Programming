appname := OMicroN

# ==========================================
# 1. 路徑設定 (關鍵修改區)
# ==========================================
# [CUDA] 指向您在家目錄安裝的 CUDA 11.4
PATH_CUDA_LOCAL := $(HOME)/cuda-11.4

PATH_HDF5 := $(HOME)/OMicroN/hdf5_local
PATH_TIFF := $(HOME)/OMicroN/libtiff_local

# ==========================================
# 2. 工具設定
# ==========================================
CXX := g++
# [CUDA] 強制使用私人安裝的 nvcc
NVCC := $(PATH_CUDA_LOCAL)/bin/nvcc

# ==========================================
# 3. 搜尋路徑設定
# ==========================================
# Header 搜尋路徑 (加入 CUDA include)
INCLUDES := -Isrc -Isrc/utils -I. -I$(PATH_HDF5)/include -I$(PATH_TIFF)/include -I$(PATH_CUDA_LOCAL)/include

# Library 搜尋路徑 (加入 CUDA lib64)
LIBS_DIRS := -L$(PATH_HDF5)/lib -L$(PATH_TIFF)/lib -L$(PATH_CUDA_LOCAL)/lib64
LIBS_NAMES := -lm -ltiff -lpthread -lhdf5_cpp -lhdf5 -ldl -lgomp

# ==========================================
# 4. 編譯參數
# ==========================================

# C++ 編譯參數
# 移除 -march=native 以避免相容性問題，改用較通用的設定
CXXFLAGS := -std=c++17 -ffast-math -Wno-unknown-pragmas -fPIC -fmax-include-depth=300 -fopenmp

# CUDA 編譯參數
# -arch=sm_61: 指定 Pascal 架構 (配合學校舊機器)
# -allow-unsupported-compiler: 忽略 GCC 版本過新的錯誤 (關鍵!)
CUDAFLAGS := -O3 -arch=sm_61 -allow-unsupported-compiler -Xcompiler -fPIC

# 連結參數 (給 nvcc 最終連結用)
# -cudart static: 強制靜態連結 (避免執行時找不到 .so)
LDFLAGS_NVCC := -cudart static -arch=sm_61 -allow-unsupported-compiler -Xcompiler -fopenmp

# ==========================================
# 5. 檔案設定
# ==========================================
SRC_DIR := src
UTILS_DIR := $(SRC_DIR)/utils

SOURCES := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(UTILS_DIR)/*.cpp)
SOURCES_TMP := $(patsubst $(SRC_DIR)/%,%,$(SOURCES))
SOURCES_TMP := $(patsubst $(UTILS_DIR)/%,%,$(SOURCES_TMP))

BIN_PAR := bin/OMicroN_parallel
BUILD_PAR := build_parallel

# 物件檔清單
OBJS_PAR  := $(patsubst %.cpp,$(BUILD_PAR)/%.o,$(SOURCES_TMP)) $(BUILD_PAR)/cuda_kernels.o

# ==========================================
# 6. 建置規則
# ==========================================

all: parallel

parallel: $(BIN_PAR)

# --- [連結] 使用私人的 NVCC 進行連結 ---
$(BIN_PAR): $(OBJS_PAR) | bin
	@echo "Linking with Local NVCC (Static CUDA Runtime)..."
	$(NVCC) $(LDFLAGS_NVCC) -o $@ $^ $(LIBS_DIRS) $(LIBS_NAMES)
	@echo "Build Success: $@"

# --- C++ 編譯 ---
$(BUILD_PAR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_PAR)
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ $<

$(BUILD_PAR)/%.o: $(UTILS_DIR)/%.cpp | $(BUILD_PAR)
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) -o $@ $<

# --- CUDA 編譯 (使用私人 NVCC) ---
$(BUILD_PAR)/%.o: $(SRC_DIR)/%.cu | $(BUILD_PAR)
	@mkdir -p $(dir $@)
	$(NVCC) $(CUDAFLAGS) $(INCLUDES) -c $< -o $@

# --- 資料夾 ---
bin $(BUILD_PAR):
	mkdir -p $@

clean:
	rm -rf bin build_parallel

.PHONY: all parallel clean