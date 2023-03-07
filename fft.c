#define _USE_MATH_DEFINES

#include <stdio.h>
#include <math.h>
#include <complex.h>

void print(const double complex* data, int size) {
    for (int i = 0; i < size; i++) {
        printf("(%g, %g) ", creal(data[i]), cimag(data[i]));
    }
    printf("\n");
}

#define N 8

int main(int argc, char* argv[]) {
    double complex data[N] = {1+I,2,3+2*I,4,5,6+I,7,8};

    print(data, N);

    double complex even[N/2];
    double complex odd[N/2];

    double complex kSum[N];

    double sumRealEven, sumImagEven, sumRealOdd, sumImagOdd; // temporary sums

    for (int k=0; k < N/2; k++) {
        printf("In loop k=%d\n", k);
        sumRealEven = 0.0;
        sumImagEven = 0.0;
        sumRealOdd = 0.0;
        sumImagOdd = 0.0;

        for (int i=0; i <= (N/2 - 1); i++) {
            printf("\tIn loop i1=%d\n", i);
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
            printf("\tIn loop i2=%d\n", i);
            sumRealEven += creal(even[i]);
            sumImagEven += cimag(even[i]);
            sumRealOdd += creal(odd[i]);
            sumImagOdd += cimag(odd[i]);
        }

        kSum[k] = (sumRealEven + sumRealOdd) + (sumImagEven + sumImagOdd)*I;
        kSum[k+N/2] = (sumRealEven - sumRealOdd) + (sumImagEven - sumImagOdd)*I;
    }
    
    print(kSum, N);

    return 0;
}