#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#include "mpi.h"

#include <stdio.h>
#include <unistd.h>

#define N_RUNS 1

void print(const double complex* data, int size) {
    for (int i = 0; i < size; i++) {
        printf("(%lf, %lf) ", creal(data[i]), cimag(data[i]));
    }
    printf("\n");
}

/*
int showInterface(const char* filename) {
    int choice = 0;
    while(!choice) {
        printf("###   SRiR - FFT   ###\n\n");
        printf("Select option:\n");

        printf("1. Read from file\n");
        printf("9. Exit\n\n");

        printf("Choice: \n");
        scanf("%d", &choice);
        // clrscr();
    }

    switch(choice) {
        case 1:
            printf("Input file path:\n");
            scanf("%s", filename);
            return 1;
        case 9:
            return 0;
        default:
            printf("No such option!\n");
            return -1;
    }
}
*/

int parse_file(const char* fileName, double complex** data) {
    FILE* file = fopen(fileName, "r");

    if (file == NULL)
        return -1;

    // read file size
    int N;
    if (fscanf(file, "%d", &N) != 1)
        return -1;

    // increase size to power of 2
    int k=2;
    while (k < N) {
        k *= 2;
    }
    N = k;

    *data = malloc(N * sizeof(double complex));
    if (data == NULL)
        return -1;

    int i = 0;
    double real, imag;
    while (i < N && fscanf(file, "%lf %lf", &real, &imag) == 2) {
        (*data)[i] = real + imag*I;
        i++;
    }

    // pad with zeros
    while (i < N) {
        (*data)[i] = 0;
        i++;
    }

    fclose(file);
    return N;
}

void write_to_file(const char* fileName, double complex* data, int N) {
    FILE* file = fopen(fileName, "w");

    if (file == NULL)
        exit(EXIT_FAILURE);

    fprintf(file, "%d\n", N);

    for (int i=0; i < N; i++) {
        fprintf(file, "%lf %lf\n", creal(data[i]), cimag(data[i]));
    }

    fclose(file);
}

void bitReversePermutation(double complex* data, int N) {
        for (int i = 0; i < N; i++) { // i - obecny indeks
        int bit_reverse_i = 0; // obecny indeks po odwróceniu kolejności bitów
        for (int j = 0; j < log2(N); j++) {
            bit_reverse_i <<= 1;
            bit_reverse_i |= (i >> j) & 1; 
        }
        if (bit_reverse_i < i) { // zamieniamy elementy zgodnie z nowymi indeksami
            double complex tmp = data[i];
            data[i] = data[bit_reverse_i];
            data[bit_reverse_i] = tmp;
        }
    }
}

void serialFFT(double complex* data, int N, bool skipPermutation, int rank) {
    if (!skipPermutation)
        bitReversePermutation(data, N);

    for (int s = 1; s <= log2(N); s++) {
        //printf("For s = %d (rank = %d)\n", s, rank);
        int m = pow(2, s);
        for (int i = 0; i < N; i += m) {
            for (int j = i; j < i + m/2; j++) {
                //printf("\tcalculating indexes: %d and %d\n", j, j+m/2);
                double complex temp = data[j];
                double complex twiddle= cexp(-2 * M_PI * (j - i) / m * I) * data[j + m/2];
                data[j] = temp + twiddle;
                data[j + m/2] = temp - twiddle;
            }
        }
    }
}

void parallelFFT(complex double* data, int N, int rank, int numProcs)
{
    if (rank == 0) // root dokonuje permutacji
        bitReversePermutation(data, N);

    int n_local = N / numProcs;
    complex double* local_data = malloc(n_local * sizeof(double complex));
    MPI_Barrier(MPI_COMM_WORLD);

    // Rozdziela dane pomiędzy procesy, ex. P1 dostanie indeksy 0-3, P2 dostanie indeksy 4-7
    MPI_Scatter(data, n_local, MPI_C_DOUBLE_COMPLEX, local_data, n_local, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Każdy proces dokonuje permutacji i liczy FFT na lokalnych danych (faza 1 i 2)
    serialFFT(local_data, n_local, true, rank);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(local_data, n_local, MPI_C_DOUBLE_COMPLEX, data, n_local, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    free(local_data);

    int Isend=0, Irecv=0;
    int c, recv_count;
    int log_n = log2(N);
    int log_p = log2(numProcs);
    for (int s = log_n - log_p + 1; s <= log_n; s++) {
        MPI_Request reqs[1000000];
        //int reqs = 0;
        int m = pow(2, s);
        c = 0;
        for (int i = 0; i < N; i += m) {
            recv_count = 0;
            for (int j = i+rank; j < i + m/2; j+=numProcs) {
                //printf("For s = %d (rank = %d) calculating indexes: %d and %d\n", s, rank, j, j+m/2);
                double complex temp = data[j];
                double complex twiddle= cexp(-2 * M_PI * (j - i) / m * I) * data[j + m/2];
                data[j] = temp + twiddle;
                data[j + m/2] = temp - twiddle;
                for (int proc=0; proc < numProcs; proc++) {
                    if (proc != rank) {
                        //printf("Sending data[%d] to proc=%d with tag=%d (from rank=%d)\n", j, proc, j, rank);
                        MPI_Isend(&data[j], 1, MPI_C_DOUBLE_COMPLEX, proc, j, MPI_COMM_WORLD, &reqs[c++]);
                        //printf("Sending data[%d] to proc=%d with tag=%d (from rank=%d)\n", j+m/2, proc, j+m/2, rank);
                        MPI_Isend(&data[j + m/2], 1, MPI_C_DOUBLE_COMPLEX, proc, j + m/2, MPI_COMM_WORLD, &reqs[c++]);
                        //printf("Receiving data[%d] from proc=%d with tag=%d (from rank=%d)\n", i+proc + recv_count*numProcs, proc, i+proc + recv_count*numProcs, rank);
                        MPI_Irecv(&data[i + proc + recv_count*numProcs], 1, MPI_C_DOUBLE_COMPLEX, proc, i + proc + recv_count*numProcs, MPI_COMM_WORLD, &reqs[c++]);
                        //printf("Receiving data[%d] from proc=%d with tag=%d (from rank=%d)\n", i+proc + recv_count*numProcs+m/2, proc, i+proc + recv_count*numProcs+m/2, rank);
                        MPI_Irecv(&data[i + proc + recv_count*numProcs + m/2], 1, MPI_C_DOUBLE_COMPLEX, proc, i + proc + recv_count*numProcs + m/2, MPI_COMM_WORLD, &reqs[c++]);
                    }
                }
                recv_count++;
            }
        }
        MPI_Waitall(c, reqs, MPI_STATUSES_IGNORE);
    }
}

int main(int argc, char* argv[]) {
    //char filename[200];
    //int code = showInterface(filename);
    //if(code != 1) return 0;

    int N;
    double complex* data;
    int rank, numProcs, nameLen;
    char   processorName[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Get_processor_name(processorName, &nameLen);

    clock_t start, end;
    double t1, t2;
    double cpu_time_used;
    
    // root process reads data from file and broadcasts
    if (rank == 0) {
        N = parse_file("data/data.txt", &data);
        printf("N = %d\n", N);
        if (N == -1) {
            printf("parse_file: something went wrong\n");
            exit(EXIT_FAILURE);
        }
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) data = malloc(N * sizeof(double complex));

    // double complex* data_copy;
    // if (rank == 0) {
    //     data_copy = malloc(N * sizeof(double complex));
    //     memcpy(data_copy, data, N * sizeof(double complex));
    // }
    
    if (rank == 0) {
        double total_time = 0.0;
        for (int i=0; i < N_RUNS; i++) {
            start = clock();
            serialFFT(data, N, false, -1);
            write_to_file("output/outputSerial.txt", data, N);
            end = clock();
            cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
            total_time += cpu_time_used;
            //memcpy(data_copy, data, N * sizeof(double complex));
            //printf("serialFFT() time taken: %lf\n", cpu_time_used);
        }
        printf("serialFFT() after %d runs: total_time: %lf, average time: %lf\n", N_RUNS, total_time, total_time/N_RUNS);
    }

    // double total_time = 0.0;
    // for (int i=0; i < N_RUNS; i++) {
    //     t1 = MPI_Wtime();
    //     parallelFFT(data, N, rank, numProcs);
    //     t2 = MPI_Wtime();
    //     if (rank == 0) {
    //         total_time += (t2-t1);
    //         printf("parallelFFT() time taken: %f\n", t2-t1);
    //         write_to_file("output/outputParallel.txt", data, N);
    //     }
    // }
    // if (rank == 0)
    //     printf("parallelFFT() after %d runs: total_time: %lf, average time: %lf\n", N_RUNS, total_time, total_time/N_RUNS);


    free(data);

    if (rank == 0) {
        printf("BEFORE FINALIZE(FINALIZE DOESN'T WORK - USE CTRL+C)\n");
    }
    MPI_Finalize();

    return 0;
}