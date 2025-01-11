#include <iostream>
#include <vector>
#include <omp.h> // OpenMP header

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
    double start_seq = omp_get_wtime();
    for (int i = 0; i < size; i++) {
        sequential_sum += arr[i];
    }
    double end_seq = omp_get_wtime();
    std::cout << "Sequential sum: " << sequential_sum
        << " computed in " << (end_seq - start_seq) << " seconds." << std::endl;

    // Parallel computation with OpenMP
    double start_parallel = omp_get_wtime();
#pragma omp parallel for reduction(+:parallel_sum)
    for (int i = 0; i < size; i++) {
        parallel_sum += arr[i];
    }
    double end_parallel = omp_get_wtime();

    // Display results
    std::cout << "Parallel sum: " << parallel_sum
        << " computed in " << (end_parallel - start_parallel) << " seconds." << std::endl;

    // Verify results
    if (parallel_sum == sequential_sum) {
        std::cout << "OpenMP test successful! Results match." << std::endl;
    }
    else {
        std::cout << "Error: Results do not match." << std::endl;
    }

    return 0;
}
