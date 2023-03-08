#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

// TODO: gui, padding data with zero, loading from/writing to file

void print(const double complex* data, int size) {
    for (int i = 0; i < size; i++) {
        printf("(%lf, %lf) ", creal(data[i]), cimag(data[i]));
    }
    printf("\n");
}

int parse_file(const char* fileName, double complex** data) {
    FILE* file = fopen(fileName, "r");

    if (file == NULL)
        return -1;

    // read file size
    int N;
    if (fscanf(file, "%d", &N) != 1)
        return -1;

    *data = malloc(N * sizeof(double complex));
    if (data == NULL)
        return -1;

    int i = 0;
    double real, imag;
    while (i < N && fscanf(file, "%lf %lf", &real, &imag) == 2) {
        (*data)[i] = real + imag*I;
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

int main() {

    double complex* data = NULL;
    int N = parse_file("data/data.txt", &data);
    if (N == -1) {
        printf("parse_file: something went wrong\n");
        exit(EXIT_FAILURE);
    }

    print(data, N);
    
    double complex even[N/2];
    double complex odd[N/2];

    double complex kSum[N];

    double sumRealEven, sumImagEven, sumRealOdd, sumImagOdd; // temporary sums

    for (int k=0; k < N/2; k++) {
        sumRealEven = 0.0;
        sumImagEven = 0.0;
        sumRealOdd = 0.0;
        sumImagOdd = 0.0;

        for (int i=0; i <= (N/2 - 1); i++) {
            // even
            double complex c1Even = data[2*i];
            double factorEven = ((2*M_PI) * ((2*i)*k)) / N;
            double complex c2Even = (cos(factorEven) - (sin(factorEven)*I));
            even[i] = c1Even * c2Even;

            // odd
            double complex c1Odd = data[2*i+1];
            double factorOdd = ((2*M_PI) * ((2*i+1)*k)) / N;
            double complex c2Odd = (cos(factorOdd) - (sin(factorOdd)*I));
            odd[i] = c1Odd * c2Odd;
        }

        for (int i=0; i < N/2; i++) {
            sumRealEven += creal(even[i]);
            sumImagEven += cimag(even[i]);
            sumRealOdd += creal(odd[i]);
            sumImagOdd += cimag(odd[i]);
        }

        kSum[k] = (sumRealEven + sumRealOdd) + (sumImagEven + sumImagOdd)*I;
        kSum[k+N/2] = (sumRealEven - sumRealOdd) + (sumImagEven - sumImagOdd)*I;
    }
    
    print(kSum, N);
    write_to_file("output/output.txt", kSum, N);

    free(data);
    return 0;
}