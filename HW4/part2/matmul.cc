#include <mpi.h>
#include <algorithm>

static int world_rank;
static int world_size;

static int row;
static int *send_count_A = nullptr;
static int *disP_A = nullptr;
static int *recv_count_C = nullptr;
static int *disP_C = nullptr;

static bool use_tiling = false;
const int tile_size = 64;

void construct_matrices(
    int n, int m, int l, const int *a_mat, const int *b_mat, int **a_mat_ptr, int **b_mat_ptr)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank == 0) {
        if ((m % 32 == 0) || (l % 32 == 0)) {
            use_tiling = true;
        }
        else {
            use_tiling = false;
        }
    }
    MPI_Bcast(&use_tiling, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    *b_mat_ptr = new int[m*l];
    if (world_rank == 0) {
        for (int i = 0; i < m*l; ++i) {
            (*b_mat_ptr)[i] = b_mat[i];
        }
    }
    MPI_Bcast(*b_mat_ptr, m*l, MPI_INT, 0, MPI_COMM_WORLD);

    int row_per_rank = n / world_size;
    int remainder = n % world_size;
    if (world_rank < remainder) {
        row = row_per_rank + 1;
    }
    else {
        row = row_per_rank;
    }
    *a_mat_ptr = new int[row*m];

    if (world_rank == 0) {
        send_count_A = new int[world_size];
        disP_A = new int[world_size];
        recv_count_C = new int[world_size];
        disP_C = new int[world_size];

        int current_disP_A = 0;
        int current_disP_C = 0;
        for (int i = 0; i < world_size; ++i) {
            int current_row_per_rank;
            if (i < remainder) {
                current_row_per_rank = row_per_rank + 1;
            }
            else {
                current_row_per_rank = row_per_rank;
            }
            
            send_count_A[i] = current_row_per_rank * m;
            disP_A[i] = current_disP_A;
            current_disP_A += send_count_A[i];

            recv_count_C[i] = current_row_per_rank * l;
            disP_C[i] = current_disP_C;
            current_disP_C += recv_count_C[i];
        }
    }
    MPI_Scatterv(a_mat, send_count_A, disP_A, MPI_INT, *a_mat_ptr, row*m, MPI_INT, 0, MPI_COMM_WORLD);
}

void matrix_multiply(
    const int n, const int m, const int l, const int *a_mat, const int *b_mat, int *out_mat)
{
    int *local_mat = new int[row*l];

    if (use_tiling) {
        for (int i = 0; i < row*l; ++i) {
            local_mat[i] = 0;
        }
        for (int ii = 0; ii < row; ii += tile_size) {
            for (int jj = 0; jj < l; jj += tile_size) {
                for (int kk = 0; kk < m; kk += tile_size) {
                    int i_end = std::min(ii + tile_size, row);
                    int j_end = std::min(jj + tile_size, l);
                    int k_end = std::min(kk + tile_size, m);

                    for (int j = jj; j < j_end; ++j) {
                        for (int i = ii; i < i_end; ++i) {
                            int partial_sum = 0;
                            for (int k = kk; k < k_end; ++k) {
                                partial_sum += a_mat[i*m + k] * b_mat[j*m + k];
                            }
                            local_mat[i*l + j] += partial_sum;
                        }
                    }
                }
            }
        }
    }
    else {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < l; ++j) {
                int sum = 0;
                for (int k = 0; k < m; ++k) {
                    sum += a_mat[i*m+k] * b_mat[j*m+k];
                }
                local_mat[i*l+j] = sum;
            }
        }
    }
    MPI_Gatherv(local_mat, row*l, MPI_INT, out_mat, recv_count_C, disP_C, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] local_mat;
}

void destruct_matrices(int *a_mat, int *b_mat)
{
    delete[] a_mat;
    delete[] b_mat;

    if (world_rank == 0) {
        delete[] send_count_A;
        delete[] disP_A;
        delete[] recv_count_C;
        delete[] disP_C;
    }
}
