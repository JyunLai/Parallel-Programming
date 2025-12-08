#include <array>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <chrono>

struct WorkerArgs
{
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int *output;
    int threadId;
    int numThreads;
};

extern void mandelbrot_serial(float x0,
                              float y0,
                              float x1,
                              float y1,
                              int width,
                              int height,
                              int start_row,
                              int num_rows,
                              int max_iterations,
                              int *output);

/* // Q1 & Q2
void worker_thread_start(WorkerArgs *const args)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    int rows_per_thread = args->height / args->numThreads;
    int remainder = args->height % args->numThreads;

    int start_row;
    int num_rows;

    if (args->threadId < remainder) {
        num_rows = rows_per_thread + 1;
        start_row = args->threadId * num_rows;
    } else {
        num_rows = rows_per_thread;
        start_row = remainder * (rows_per_thread + 1) + (args->threadId - remainder) * rows_per_thread;
    }

    mandelbrot_serial(args->x0, args->y0, args->x1, args->y1,
                      args->width, args->height,
                      start_row, num_rows,
                      args->maxIterations, args->output);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    
    printf("Thread %d took %.4f ms to complete.\n", args->threadId, duration.count());
}
*/
/*
// Q3 & Q4
void worker_thread_start(WorkerArgs *const args)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    for (unsigned int row = args->threadId; row < args->height; row += args->numThreads){
        mandelbrot_serial(args->x0, args->y0, args->x1, args->y1,
                          args->width, args->height,
                          row,
                          1,
                          args->maxIterations, args->output);
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    printf("Thread %d (Interleaved) took %.4f ms to complete.\n", args->threadId, duration.count());
    
}
*/
// 最終優化版本 (手動內聯 + 微觀優化)
void worker_thread_start(WorkerArgs *const args)
{
    // 1. 將常數計算提到最外層，只計算一次
    float dx = (args->x1 - args->x0) / args->width;
    float dy = (args->y1 - args->y0) / args->height;

    // 2. 外層迴圈依然是高效的交錯式分配
    for (unsigned int y = args->threadId; y < args->height; y += args->numThreads) {
        
        // --- 手動內聯 mandelbrot_serial 的邏輯 ---
        // 預先計算好當前行的 c_im 和在 output 陣列中的起始位置
        float c_im = args->y0 + y * dy;
        int row_start_index = y * args->width;

        // 3. 內層迴圈處理該行的每一個像素
        for (unsigned int x = 0; x < args->width; ++x) {
            
            float c_re = args->x0 + x * dx;

            // --- 手動內聯 mandel 函式的邏輯 ---
            float z_re = c_re;
            float z_im = c_im;
            int iter;

            // 4. 這是計算最密集的「熱迴圈」
            for (iter = 0; iter < args->maxIterations; ++iter) {
                
                // 💡 *** 關鍵微觀優化 ***
                // 預先計算平方值，避免在下一行重複計算
                float z_re_sq = z_re * z_re;
                float z_im_sq = z_im * z_im;

                // 判斷是否逃逸
                if (z_re_sq + z_im_sq > 4.f) {
                    break;
                }

                // 計算下一次迭代的值 (使用已算好的平方值)
                float new_re = z_re_sq - z_im_sq;
                float new_im = 2.f * z_re * z_im;
                z_re = c_re + new_re;
                z_im = c_im + new_im;
            }
            
            // 5. 直接將結果寫入 output 陣列
            args->output[row_start_index + x] = iter;
        }
    }
}

//
// mandelbrot_thread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
void mandelbrot_thread(int num_threads,
                       float x0,
                       float y0,
                       float x1,
                       float y1,
                       int width,
                       int height,
                       int max_iterations,
                       int *output)
{
    static constexpr int max_threads = 32;

    if (num_threads > max_threads)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", max_threads);
        exit(1);
    }

    // Creates thread objects that do not yet represent a thread.
    std::array<std::thread, max_threads> workers;
    std::array<WorkerArgs, max_threads> args = {};

    for (int i = 0; i < num_threads; i++)
    {
        // TODO FOR PP STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = max_iterations;
        args[i].numThreads = num_threads;
        args[i].output = output;

        args[i].threadId = i;
    }

    // Spawn the worker threads.  Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
    for (int i = 1; i < num_threads; i++)
    {
        workers[i] = std::thread(worker_thread_start, &args[i]);
    }

    worker_thread_start(&args[0]);

    // join worker threads
    for (int i = 1; i < num_threads; i++)
    {
        workers[i].join();
    }
}
