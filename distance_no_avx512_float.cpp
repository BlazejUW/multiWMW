#include <omp.h> // OpenMP header
#include <Rcpp.h>
#include <cmath>
#include <memory>
#include <iostream>
#include <cstdlib> // for posix_memalign and free

using namespace Rcpp;

// Funkcja alokująca wyrównaną pamięć
std::unique_ptr<float[], decltype(&std::free)> aligned_alloc_float(size_t size) {
    float* ptr;
    if (posix_memalign(reinterpret_cast<void**>(&ptr), 64, size * sizeof(float)) != 0) {
        return std::unique_ptr<float[], decltype(&std::free)>(nullptr, std::free);
    }
    return std::unique_ptr<float[], decltype(&std::free)>(ptr, std::free);
}


// Funkcja obliczająca odległość euklidesową bez użycia AVX-512 dla float
float euclidean_distance_no_avx_float(const float* a, const float* b, size_t length) {
    float sum = 0.0f;
    for (size_t i = 0; i < length; ++i) {
        float diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

// [[Rcpp::export]]
NumericMatrix dist_no_avx_parallel_float(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    auto aligned_Z = aligned_alloc_float(n * dim);

    if (!aligned_Z) {
        stop("Allocation failed");
    }

    // Kopiowanie danych do wyrównanych buforów
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            aligned_Z[i * dim + j] = static_cast<float>(Z(i, j));
        }
    }
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const float* vec_i = aligned_Z.get() + i * dim;
            const float* vec_j = aligned_Z.get() + j * dim;

            float dist = euclidean_distance_no_avx_float(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;
    }

    return DZ;
}


// [[Rcpp::export]]
NumericMatrix dist_no_avx_no_parallel_float(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    auto aligned_Z = aligned_alloc_float(n * dim);

    if (!aligned_Z) {
        stop("Allocation failed");
    }

    // Kopiowanie danych do wyrównanych buforów
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            aligned_Z[i * dim + j] = static_cast<float>(Z(i, j));
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const float* vec_i = aligned_Z.get() + i * dim;
            const float* vec_j = aligned_Z.get() + j * dim;

            float dist = euclidean_distance_no_avx_float(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;
    }

    return DZ;
}