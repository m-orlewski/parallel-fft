#include <cstdio>
#include <cstdlib>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <upcxx/upcxx.hpp>

namespace global
{
int process_count;
int slave_count;
int rank;

std::vector<float> data;
int data_size;
}

int reverseBits(int number, int range)
{
    int reversed = 0;
    for(int i = 0; i < range; i++)
    {
        reversed |= ((number >> i) & 1) << (range - 1 - i);
    }
    return reversed;
}

void parseFile(const char *path)
{
    upcxx::global_ptr<float> data = nullptr;
    upcxx::global_ptr<int> data_size = nullptr;

    if(global::rank == 0)
    {
        std::vector<float> buffer;

        std::ifstream file(path);

        if(not file.is_open())
        {
            std::cout << "Error opening file" << std::endl;
            return;
        }

        buffer.push_back(0);

        float val{};

        while(file >> val)
            buffer.push_back(val);
        
        data = upcxx::new_array<float>(buffer.size());
        data_size = upcxx::new_<int>();

        for(int i = 0; i < buffer.size(); i++)
            upcxx::rput(buffer[i], data + i).wait();

        upcxx::rput(static_cast<int>(buffer.size()), data_size).wait();
    }

    data = upcxx::broadcast(data, 0).wait();
    data_size = upcxx::broadcast(data_size, 0).wait();

    global::data_size = upcxx::rget(data_size).wait();

    for(int i = 0; i < global::data_size; i++)
        global::data.push_back(upcxx::rget(data + i).wait());
}

void printResults(const float *tab_re, const float *tab_im, double time)
{
    if(0 == global::rank)
    {
        for(int i = 1; i < global::data_size; i++)
            std::cout << tab_re[i] << "\t" << tab_im[i] << " i" << std::endl;

        std::cout << "Run time: " << time << " ms" << std::endl;
    }
}

int main()
{
    upcxx::init();
    
    global::process_count = upcxx::rank_n();
    global::slave_count = global::process_count - 1;
    global::rank = upcxx::rank_me();

    parseFile("data/data.txt");

    upcxx::global_ptr<float> tab_re = nullptr;
    upcxx::global_ptr<float> tab_im = nullptr;

    if(global::rank == 0)
    {
        tab_re = upcxx::new_array<float>(global::data_size);
        tab_im = upcxx::new_array<float>(global::data_size);

        const int max_bit_width = std::log2f(global::data_size);
        for(int i = 0; i < global::data_size; i++)
        {
            upcxx::rput(global::data[reverseBits(i - 1, max_bit_width) + 1], tab_re + i).wait();
            upcxx::rput(0.0f, tab_im + i).wait();
        }
    }

    tab_re = upcxx::broadcast(tab_re, 0).wait();
    tab_im = upcxx::broadcast(tab_im, 0).wait();

    const int size_local = global::data_size / global::slave_count;

    auto start_time = std::chrono::high_resolution_clock::now();
    for(int div = 1, key = std::log2f(global::data_size - 1); key > 0; key--, div *= 2)
    {
        float buffer_re[size_local]{};
        float buffer_im[size_local]{};
        if(global::rank != 0)
        {
            for(int i = 0; i < size_local; i++)
            {
                const auto offset = (global::rank - 1) * size_local + i + 1;
                const auto is_even = ((offset + div - 1) / div) % 2;
                const auto is_odd = 1 - is_even;
                const auto butterfly_index = M_PI * ((offset - 1) % (div * 2)) / div;

                float val_re_odd = upcxx::rget(tab_re + offset - (div * is_odd)).wait();
                float val_re_even = upcxx::rget(tab_re + offset + (div * is_even)).wait();
                float val_im_odd = upcxx::rget(tab_im + offset - (div * is_odd)).wait();
                float val_im_even = upcxx::rget(tab_im + offset + (div * is_even)).wait();

                buffer_re[i] = val_re_odd + (std::cos(butterfly_index) * val_re_even)
                               + (std::sin(butterfly_index) * val_im_even);

                buffer_im[i] = val_im_odd + (std::cos(butterfly_index) * val_im_even)
                              - (std::sin(butterfly_index) * val_re_even);
            }
        }

        upcxx::barrier();

        if(global::rank != 0)
        {
            for(int i = 0; i < size_local; i++)
            {
                const auto offset = (global::rank - 1) * size_local + i + 1;
                upcxx::rput(buffer_re[i], tab_re + offset).wait();
                upcxx::rput(buffer_im[i], tab_im + offset).wait();
            }
        }

        upcxx::barrier();
    }
    auto end_time = std::chrono::high_resolution_clock::now();

    if(global::rank == 0)
    {
        double time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        float *tab_re_arr = tab_re.local();
        float *tab_im_arr = tab_im.local();
        printResults(tab_re_arr, tab_im_arr, time);
    }

    upcxx::finalize();
    return 0;
}