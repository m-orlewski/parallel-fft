#include <cstdio>
#include <cstdlib>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <complex>
#include <upcxx/upcxx.hpp>

const double pi = std::acos(-1.0);

void parse_file(std::string fileName, std::vector<std::complex<double>>& input_data, int& N) {
    FILE* file = fopen(fileName.c_str(), "r");

    if (file == NULL) {
        std::cerr << "Nie udało się otworzyć pliku do wczytania danych" << std::endl;
        return;
    }

    // Zwiększa N do najbliższej potęgi 2
    int k=2;
    while (k < N) {
        k *= 2;
    }
    if (k > N) {
        std::cout << "Ilość wczytanych danych zwiększono z " << N << " do " << k << std::endl;
        N = k;
    }

    // Wczytaj dane z pliku
    int i = 0;
    double real, imag;
    while (i < N && fscanf(file, "%lf %lf", &real, &imag) == 2) {
        input_data.push_back({real, imag});
        i++;
    }

    // Uzupełnij zerami
    while (i < N) {
        input_data.push_back({0.0, 0.0});
        i++;
    }

    fclose(file);
}

void write_to_file(std::string fileName, upcxx::global_ptr<std::complex<double>>& data, int N) {
    FILE* file = fopen(fileName.c_str(), "w");

    if (file == NULL) {
        std::cerr << "Nie udało się otworzyć pliku do zapisu danych" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i=0; i < N; i++) {
        fprintf(file, "%lf %lf\n", data.local()[i].real(), data.local()[i].imag());
    }

    fclose(file);
}

void bitReversePermutation(std::vector<std::complex<double>>& data, int N) {
    for (int i = 0; i < N; i++) { // i - obecny indeks
        int bit_reverse_i = 0; // obecny indeks po odwróceniu kolejności bitów
        for (int j = 0; j < log2(N); j++) {
            bit_reverse_i <<= 1;
            bit_reverse_i |= (i >> j) & 1; 
        }
        if (bit_reverse_i < i) { // zamieniamy elementy zgodnie z nowymi indeksami
            auto tmp = data[i];
            data[i] = data[bit_reverse_i];
            data[bit_reverse_i] = tmp;
        }
    }
}

void serialFFT(upcxx::global_ptr<std::complex<double>>& data, int N, int rank) {
    for (int s = 1; s <= log2(N); s++) {
        int m = pow(2, s);
        for (int i = 0; i < N; i += m) {
            for (int j = i; j < i + m/2; j++) {
                std::complex<double> p(0, -2 * pi * (j - i) / m);

                auto temp = upcxx::rget(data + j).wait();
                auto twiddle = std::exp(p) * upcxx::rget(data + j + m/2).wait();

                upcxx::rput(temp + twiddle, data + j).wait();
                upcxx::rput(temp - twiddle, data + j + m/2).wait();
            }
        }
    }
}

int main(int argc, char* argv[])
{
    upcxx::init();

    int numProcs = upcxx::rank_n();
    int rank = upcxx::rank_me();
    int N;
    std::string filename;

    if (argc < 3) {
        if (rank == 0)
            std::cerr << "Podaj argumenty wywołania programu: {rozmiar FFT} {ścieżka do pliku z danymi}" << std::endl;
        upcxx::finalize();
        return -1;
    }
    else {
        N = atoi(argv[1]);
        filename = argv[2];

        if (N % numProcs) {
            std::cerr << "Rozmiar FFT musi być podzielny przez numProcs" << std::endl;
        }
    }

    // Root wczytuje dane i broadcastuje je do pozostałych procesów
    std::vector<std::complex<double>> input_data;
    upcxx::global_ptr<std::complex<double>> data = nullptr;
    if (rank == 0) {
        std::cout << "Wczytywanie " << N << " danych z pliku " << filename << std::endl;
        parse_file(filename, input_data, N);

        bitReversePermutation(input_data, N);

        data = upcxx::new_array<std::complex<double>>(N);
        for (int i=0; i < N; i++) {
            upcxx::rput(input_data[i], data+i).wait();
        }
    }
    
    data = upcxx::broadcast(data, 0).wait();
    upcxx::barrier();

    int n_local = N/numProcs;
    upcxx::global_ptr<std::complex<double>> local_data = data + rank * n_local;

    // Każdy proces wykonuje log(N) - log(numProcs) iteracji
    serialFFT(local_data, n_local, rank);
    upcxx::barrier();

    // Root dokańcza obliczenia (log(numProcs) iteracji)
    if (rank == 0 && data.is_local()) {
        int log_n = log2(N);
        int log_p = log2(numProcs);
        for (int s = log_n - log_p + 1; s <= log_n; s++) {
            int m = pow(2, s);
            for (int i = 0; i < N; i += m) {
                for (int j = i; j < i + m/2; j++) {
                    std::complex<double> p(0, -2 * pi * (j - i) / m);

                    auto temp = data.local()[j];
                    auto twiddle = std::exp(p) * data.local()[j + m/2];

                    data.local()[j] = temp + twiddle;
                    data.local()[j + m/2] = temp - twiddle;
                }
            }
        }
    }

    if (rank == 0) {
        write_to_file("output/output.txt", data, N);
        std::cout << "Zapisano wyniki do pliku output/output.txt" << std::endl;
    }


    upcxx::finalize();
    return 0;
}