#include <iostream>
#include <vector>
#include <chrono> // For high-resolution timing
#include <omp.h>  // For OpenMP

int main() {
    // Array size
    const int size = 100000000; // 100 million elements
    std::vector<int> arr(size, 1); // Initialize an array of 1s

    // Variables for sum
    long long parallel_sum = 0;
    long long sequential_sum = 0;

    // Set number of threads
    int num_threads = 4; // You can modify this to test different thread counts
    omp_set_num_threads(num_threads);

    // Sequential computation (for comparison)
    auto start_seq = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < size; i++) {
        sequential_sum += arr[i];
    }
    auto end_seq = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seq = end_seq - start_seq;
    std::cout << "Sequential sum: " << sequential_sum
        << " computed in " << elapsed_seq.count() << " seconds." << std::endl;

    // Parallel computation with OpenMP
    auto start_parallel = std::chrono::high_resolution_clock::now();
#pragma omp parallel for reduction(+:parallel_sum)
    for (int i = 0; i < size; i++) {
        parallel_sum += arr[i];
    }
    auto end_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel = end_parallel - start_parallel;

    // Display results
    std::cout << "Parallel sum: " << parallel_sum
        << " computed in " << elapsed_parallel.count() << " seconds." << std::endl;

    // Verify results
    if (parallel_sum == sequential_sum) {
        std::cout << "OpenMP test successful! Results match." << std::endl;
    }
    else {
        std::cout << "Error: Results do not match." << std::endl;
    }

    return 0;
}
