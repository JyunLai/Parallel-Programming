/* src/cuda_kernels.cu */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>
#include <stdio.h>
#include <vector>

// 定義常數
#ifndef M_PI
#define M_PI 3.14159265359f
#endif

#define MAX_NEIGHBORS 26
#define N_SYM 24 // 立方晶系有 24 個對稱操作

// =========================================================
// 1. Device 端的基本數學結構與函式 (取代 Eigen)
// =========================================================

// 簡單的四元數結構 (w, x, y, z)
struct float4_quat {
    float w, x, y, z;
};

// 簡單的向量結構
struct float3_vec {
    float x, y, z;
};

__device__ inline float d2r(float d) { return d * M_PI / 180.0f; }

// 四元數乘法: q = q1 * q2
__device__ float4_quat quat_mult(float4_quat q1, float4_quat q2) {
    float4_quat res;
    res.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    res.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    res.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    res.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return res;
}

// 四元數共軛
__device__ float4_quat quat_conj(float4_quat q) {
    return {q.w, -q.x, -q.y, -q.z};
}

// 四元數正規化
__device__ float4_quat quat_normalize(float4_quat q) {
    float norm = sqrtf(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
    if (norm > 1e-9f) {
        float inv = 1.0f / norm;
        return {q.w * inv, q.x * inv, q.y * inv, q.z * inv};
    }
    return {1.0f, 0.0f, 0.0f, 0.0f}; // Default identity
}

// Euler Angles (phi1, Phi, phi2) 轉 Quaternion
// 對應 orientation.cpp 中的 ea2q
__device__ float4_quat ea2q_device(float phi1, float Phi, float phi2) {
    float p1 = d2r(phi1) * 0.5f;
    float P  = d2r(Phi)  * 0.5f;
    float p2 = d2r(phi2) * 0.5f;

    float c1 = cosf(p1), s1 = sinf(p1);
    float cP = cosf(P),  sP = sinf(P);
    float c2 = cosf(p2), s2 = sinf(p2);

    // Z1 * X * Z2 順序 (Bunge Convention)
    // q = [c1 cP c2 - s1 cP s2, c1 sP c2 + s1 sP s2, c1 sP s2 - s1 sP c2, c1 cP s2 + s1 cP c2]
    // 注意：Eigen 的 Quaternion 建構順序通常是 W, X, Y, Z，但內部儲存順序可能不同。
    // 這裡手刻標準轉換公式：
    
    float4_quat q;
    q.w = c1*cP*c2 - s1*cP*s2;
    q.x = c1*sP*c2 + s1*sP*s2;
    q.y = s1*sP*c2 - c1*sP*s2;
    q.z = c1*cP*s2 + s1*cP*c2;
    
    return quat_normalize(q);
}

// =========================================================
// 2. 立方晶系對稱操作 (Hardcoded in Constant Memory)
// =========================================================
// 這是立方晶系的 24 個對稱四元數，從 orientation.cpp 的邏輯推導而來
// 為了效能，我們直接寫死在 GPU 程式碼中
__constant__ float4_quat d_SYM[24] = {
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1},
    {0.5, 0.5, 0.5, 0.5}, {0.5, 0.5, 0.5, -0.5}, {0.5, 0.5, -0.5, 0.5}, {0.5, 0.5, -0.5, -0.5},
    {0.5, -0.5, 0.5, 0.5}, {0.5, -0.5, 0.5, -0.5}, {0.5, -0.5, -0.5, 0.5}, {0.5, -0.5, -0.5, -0.5},
    {0.70710678, 0.70710678, 0, 0}, {0.70710678, -0.70710678, 0, 0}, 
    {0.70710678, 0, 0.70710678, 0}, {0.70710678, 0, -0.70710678, 0},
    {0.70710678, 0, 0, 0.70710678}, {0.70710678, 0, 0, -0.70710678},
    {0, 0.70710678, 0.70710678, 0}, {0, 0.70710678, -0.70710678, 0},
    {0, 0.70710678, 0, 0.70710678}, {0, 0.70710678, 0, -0.70710678},
    {0, 0, 0.70710678, 0.70710678}, {0, 0, 0.70710678, -0.70710678}
};

// 計算兩個 Euler Angles 之間的 Misorientation
// 對應 Orientations::CalcMisorientationFromEulerAngles
__device__ float d_GetMisorientation(float ea1_0, float ea1_1, float ea1_2, 
                                     float ea2_0, float ea2_1, float ea2_2) 
{
    float4_quat q0 = ea2q_device(ea1_0, ea1_1, ea1_2);
    float4_quat q1 = ea2q_device(ea2_0, ea2_1, ea2_2);
    
    float m = fabs(q0.w*q1.w + q0.x*q1.x + q0.y*q1.y + q0.z*q1.z);
    
    // 遍歷 24 個對稱操作 (省略 i=0 因為上面已經算過 Identity)
    for (int i = 1; i < N_SYM; ++i) {
        float4_quat sym = d_SYM[i];
        // qs = q1 * SYM[i]
        float4_quat qs = quat_mult(q1, sym);
        
        float w = fabs(q0.w*qs.w + q0.x*qs.x + q0.y*qs.y + q0.z*qs.z);
        if (w > m) m = w;
    }
    
    // 防止 acos 越界
    if (m > 1.0f) m = 1.0f;
    
    return 180.0f * 2.0f * acosf(m) / M_PI;
}

// =========================================================
// 3. 物理計算輔助函式
// =========================================================

__device__ double d_GetBoundaryEnergy(double mis, double HAGB, double LowerMisCutOff, double GrainBoundaryEnergy) {
    if (mis < LowerMisCutOff) return 0.0;
    double r = mis / HAGB;
    if (r >= 1.0) return GrainBoundaryEnergy;
    return GrainBoundaryEnergy * r * (1.0 - log(r));
}

__device__ double d_GetMobility(double mis, double HAGB, double LowerMisCutOff, double MaxMobility) {
    if (mis < LowerMisCutOff) return 0.0;
    double r = mis / HAGB;
    if (r >= 1.0) return MaxMobility;
    
    double HumphreysB = 5.0;
    double HumphreysN = 4.0;
    return MaxMobility * (1.0 - exp(-HumphreysB * pow(r, HumphreysN)));
}

// 取得 3D 網格座標的 Index (假設 Periodic Boundary 為 false)
__device__ int d_GetNeighborIndex(int idx, int dx, int dy, int dz, int Nx, int Ny, int Nz) {
    int z = idx / (Nx * Ny);
    int rem = idx % (Nx * Ny);
    int y = rem / Nx;
    int x = rem % Nx;

    int nx = x + dx;
    int ny = y + dy;
    int nz = z + dz;

    // 邊界檢查 (Non-periodic)
    if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 || nz >= Nz) return -1;

    return nz * Nx * Ny + ny * Nx + nx;
}


// =========================================================
// 4. CUDA Kernel
// =========================================================

__global__ void ReorientationKernel(
    int numInterfaceCells,
    const int* __restrict__ d_indices,       // 待處理的 Cell Indices
    const int* __restrict__ d_latticeId,     // 狀態: Lattice ID
    const int* __restrict__ d_oriId,         // 狀態: Ori ID
    const float* __restrict__ d_eulerAngles, // 取向庫: [id*3, id*3+1, id*3+2]
    const float* __restrict__ d_ci,          // 狀態: CI
    const int* __restrict__ d_rx,            // 狀態: RX
    float* d_consumptionRate,                // [輸出] 速率
    int* d_neighborGrowing,                  // [輸出] 生長來源
    
    // 參數包
    int targetLatticeId,
    int Nx, int Ny, int Nz, float Dx,
    float MinCI, float LowerMisCutOff, float LowerMisForLAGB, float HAGB, 
    float GBEnergy, float Mobility, float BetaFactor
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numInterfaceCells) return;

    int index = d_indices[idx]; // 真實 Cell Index

    // 1. 檢查相 (Lattice)
    if (d_latticeId[index] != targetLatticeId) return;

    // 初始化變數
    int myOriId = d_oriId[index];
    float myEA[3] = {
        d_eulerAngles[myOriId * 3 + 0],
        d_eulerAngles[myOriId * 3 + 1],
        d_eulerAngles[myOriId * 3 + 2]
    };
    
    float fDot = 0.0f;
    int bestNeighbor = -1;
    double CurrentEnergy = 0.0;

    // 定義 26 鄰居的偏移 (簡易版，不含 GridType=1 的六角網格處理，僅作示範)
    // 您可能需要根據 settings.cpp 裡的 GetNeighbourCellOffsets 調整
    int neighborOffsets[26][3];
    int nCount = 0;
    for(int kz=-1; kz<=1; kz++) {
        for(int ky=-1; ky<=1; ky++) {
            for(int kx=-1; kx<=1; kx++) {
                if(kx==0 && ky==0 && kz==0) continue;
                if(Nz==1 && kz!=0) continue; // 2D 模式
                neighborOffsets[nCount][0]=kx; 
                neighborOffsets[nCount][1]=ky; 
                neighborOffsets[nCount][2]=kz;
                nCount++;
            }
        }
    }

    // --- Loop 1: 計算當前狀態能量 ---
    for (int n = 0; n < nCount; n++) {
        int nbrIdx = d_GetNeighborIndex(index, neighborOffsets[n][0], neighborOffsets[n][1], neighborOffsets[n][2], Nx, Ny, Nz);
        if (nbrIdx == -1) continue;

        int nbrOriId = d_oriId[nbrIdx];
        if (myOriId == nbrOriId) continue; // 相同取向略過

        // 取得鄰居角度
        float nbrEA[3] = {
            d_eulerAngles[nbrOriId * 3 + 0],
            d_eulerAngles[nbrOriId * 3 + 1],
            d_eulerAngles[nbrOriId * 3 + 2]
        };

        float mis = d_GetMisorientation(myEA[0], myEA[1], myEA[2], nbrEA[0], nbrEA[1], nbrEA[2]);
        
        if (mis < LowerMisCutOff) continue;

        // 計算能量
        double sigma = d_GetBoundaryEnergy(mis, HAGB, LowerMisCutOff, GBEnergy);
        
        // 這裡假設所有邊界權重為 1 (實際應根據 GridType 計算面積)
        double area = 1.0; 
        CurrentEnergy += area * sigma;
    }

    // --- Loop 2: 計算候選鄰居的驅動力 ---
    for (int n = 0; n < nCount; n++) {
        int nbrIdx = d_GetNeighborIndex(index, neighborOffsets[n][0], neighborOffsets[n][1], neighborOffsets[n][2], Nx, Ny, Nz);
        if (nbrIdx == -1) continue;

        // 檢查鄰居是否為合格的生長候選者 (CI, RX 狀態等)
        // (這裡簡化判斷邏輯)
        if (d_rx[nbrIdx] == 0 && d_ci[nbrIdx] < MinCI) continue;

        int candidateOriId = d_oriId[nbrIdx];
        float candEA[3] = {
            d_eulerAngles[candidateOriId * 3 + 0],
            d_eulerAngles[candidateOriId * 3 + 1],
            d_eulerAngles[candidateOriId * 3 + 2]
        };

        float misWithMe = d_GetMisorientation(myEA[0], myEA[1], myEA[2], candEA[0], candEA[1], candEA[2]);
        if (misWithMe < LowerMisForLAGB) continue; // 移動力太低

        // 計算如果變成這個鄰居，能量會變多少 (NextEnergy)
        double NextEnergy = 0.0;
        
        // Loop 2.1: 對這個候選取向，再次遍歷所有鄰居
        for (int nn = 0; nn < nCount; nn++) {
            int otherNbrIdx = d_GetNeighborIndex(index, neighborOffsets[nn][0], neighborOffsets[nn][1], neighborOffsets[nn][2], Nx, Ny, Nz);
            if (otherNbrIdx == -1) continue;
            
            int otherOriId = d_oriId[otherNbrIdx];
            if (candidateOriId == otherOriId) continue;

             float otherEA[3] = {
                d_eulerAngles[otherOriId * 3 + 0],
                d_eulerAngles[otherOriId * 3 + 1],
                d_eulerAngles[otherOriId * 3 + 2]
            };

            float misNext = d_GetMisorientation(candEA[0], candEA[1], candEA[2], otherEA[0], otherEA[1], otherEA[2]);
            
            if (misNext < LowerMisCutOff) continue;
            
            double sigma = d_GetBoundaryEnergy(misNext, HAGB, LowerMisCutOff, GBEnergy);
            double area = 1.0;
            NextEnergy += area * sigma;
        }

        double dG = CurrentEnergy - NextEnergy;
        if (dG <= 0.0) continue; // 沒有驅動力

        // 簡化：Force = dG
        double totalForce = dG; 
        
        double mob = d_GetMobility(misWithMe, HAGB, LowerMisCutOff, Mobility);
        double rate = totalForce * mob * BetaFactor;

        if (rate > fDot) {
            fDot = rate;
            bestNeighbor = nbrIdx;
        }
    }

    // 寫回結果
    d_consumptionRate[index] = fDot;
    d_neighborGrowing[index] = bestNeighbor;
}

// =========================================================
// 5. Host Wrapper
// =========================================================

extern "C" void LaunchCudaReorientation(
    int numInterfaceCells,
    const std::vector<int>& h_indices,
    int numTotalCells,
    // State Arrays
    const int* h_latticeId,
    const int* h_oriId,
    const float* h_eulerAngles, // Flattened Euler Angles [N_ORIS * 3]
    int numOris,                // 總取向數
    const float* h_ci,
    const int* h_rx,
    // Outputs
    float* h_consumptionRate_out,
    int* h_neighborGrowing_out,
    // Parameters
    int targetLatticeId,
    int Nx, int Ny, int Nz, float Dx,
    float MinCI, float LowerMisCutOff, float LowerMisForLAGB, float HAGB, 
    float GBEnergy, float Mobility, float BetaFactor
) {
    // 1. Allocate Device Memory
    int *d_indices, *d_latticeId, *d_oriId, *d_rx, *d_neighborOut;
    float *d_eulerAngles, *d_ci, *d_rateOut;

    // 錯誤檢查宏 (建議加上)
    #define CHECK(call) { cudaError_t err = call; if (err != cudaSuccess) printf("CUDA Error: %s\n", cudaGetErrorString(err)); }

    CHECK(cudaMalloc(&d_indices, numInterfaceCells * sizeof(int)));
    CHECK(cudaMalloc(&d_latticeId, numTotalCells * sizeof(int)));
    CHECK(cudaMalloc(&d_oriId, numTotalCells * sizeof(int)));
    CHECK(cudaMalloc(&d_eulerAngles, numOris * 3 * sizeof(float)));
    CHECK(cudaMalloc(&d_ci, numTotalCells * sizeof(float)));
    CHECK(cudaMalloc(&d_rx, numTotalCells * sizeof(int)));
    CHECK(cudaMalloc(&d_rateOut, numTotalCells * sizeof(float)));
    CHECK(cudaMalloc(&d_neighborOut, numTotalCells * sizeof(int)));

    // 2. Copy Host to Device
    CHECK(cudaMemcpy(d_indices, h_indices.data(), numInterfaceCells * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_latticeId, h_latticeId, numTotalCells * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_oriId, h_oriId, numTotalCells * sizeof(int), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_eulerAngles, h_eulerAngles, numOris * 3 * sizeof(float), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_ci, h_ci, numTotalCells * sizeof(float), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_rx, h_rx, numTotalCells * sizeof(int), cudaMemcpyHostToDevice));

    // 3. Launch Kernel
    int blockSize = 256;
    int gridSize = (numInterfaceCells + blockSize - 1) / blockSize;

    ReorientationKernel<<<gridSize, blockSize>>>(
        numInterfaceCells, d_indices, d_latticeId, d_oriId, d_eulerAngles,
        d_ci, d_rx, d_rateOut, d_neighborOut,
        targetLatticeId, Nx, Ny, Nz, Dx, MinCI, LowerMisCutOff, LowerMisForLAGB,
        HAGB, GBEnergy, Mobility, BetaFactor
    );

    cudaDeviceSynchronize();

    // 4. Copy Back
    // 注意：這裡只複製全部陣列回來最簡單，雖然可能有點浪費
    CHECK(cudaMemcpy(h_consumptionRate_out, d_rateOut, numTotalCells * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(h_neighborGrowing_out, d_neighborOut, numTotalCells * sizeof(int), cudaMemcpyDeviceToHost));

    // 5. Free
    cudaFree(d_indices); cudaFree(d_latticeId); cudaFree(d_oriId);
    cudaFree(d_eulerAngles); cudaFree(d_ci); cudaFree(d_rx);
    cudaFree(d_rateOut); cudaFree(d_neighborOut);
}