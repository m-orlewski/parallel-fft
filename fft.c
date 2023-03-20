#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <time.h>

#include "mpi.h"

void print(const double complex* data, int size) {
    for (int i = 0; i < size; i++) {
        printf("(%lf, %lf) ", creal(data[i]), cimag(data[i]));
    }
    printf("\n");
}

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

int logBase2(int N)    /*function to calculate the log2(.) of int numbers*/
{
    int k = N, i = 0;
    while(k) {
        k >>= 1;
        i++;
    }
    return i - 1;
}

int reverse_bits(int N, int n)    //calculating revers number
{
  int j, p = 0;
  for(j = 1; j <= logBase2(N); j++) {
    if(n & (1 << (logBase2(N) - j)))
      p |= 1 << (j - 1);
  }
  return p;
}

void reorder(complex double* f1, int N) //using the reverse order in the array
{
  complex double f2[200];
  for(int i = 0; i < N; i++)
    f2[i] = f1[reverse_bits(N, i)];
  for(int j = 0; j < N; j++)
    f1[j] = f2[j];
}

void transform(complex double* data, int N)
{
    //printf("BEFORE REORDER: \n");
    //print(data, N);
    reorder(data, N);    //first: reverse order
    //printf("AFTER REORDER: \n");
    //print(data, N);
    complex double *W;
    W = (complex double *)malloc(N / 2 * sizeof(complex double));
    W[1] = 1.0 * (cos(-2.0*M_PI/N) + sin(-2.0*M_PI/N)*I);
    //printf("W[1] = %lf, %lf\n", creal(W[1]), cimag(W[1]));
    W[0] = 1;
    for(int i = 2; i < N / 2; i++) {
        //W[i] = pow(W[1], i);
        W[i] = W[i-1]*W[1];
    }
    int n = 1;
    int a = N / 2;
    for(int j = 0; j < logBase2(N); j++) {
        for(int i = 0; i < N; i++) {
            if(!(i & n)) {
                complex double temp = data[i];
                complex double Temp = W[(i * a) % (n * a)] * data[i + n];
                data[i] = temp + Temp;
                data[i + n] = temp - Temp;
                //printf("data[i]=%lf, %lf\n", creal(data[i]), cimag(data[i]));
            }
        }
        n *= 2;
        a = a / 2;
    }
    free(W);
}

double serial_fft(double complex* data, int N) {
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    transform(data, N);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("serial_fft() time taken: %lf\n", cpu_time_used);
    
    write_to_file("output/output.txt", data, N);
    return cpu_time_used;
}

/*
void parallel_fft(double complex* data, int N) {
    Mpi_Init();

    Mpi_Finalize();
}
*/

int main(int argc, char* argv[]) {
    //char* filename = malloc(sizeof(char)*200);
    //int code = showInterface(filename);
    //if(code != 1) return 0;

    double complex* data = NULL;
    int N = parse_file("data/data.txt", &data);
    if (N == -1) {
        printf("parse_file: something went wrong\n");
        exit(EXIT_FAILURE);
    }

    //free(filename);

    //printf("%d\n", N);
    //print(data, N);
    
    double serial_time = serial_fft(data, N);
    //double parallel_time = parallel_fft(data, N);

    free(data);
    return 0;
}