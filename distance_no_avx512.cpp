#include <omp.h> // OpenMP header
#include <Rcpp.h>
#include <cmath>
#include <memory>
#include <iostream>
#include <cstdlib> // for posix_memalign and free

using namespace Rcpp;

// Funkcja alokująca wyrównaną pamięć
std::unique_ptr<double[], decltype(&std::free)> aligned_alloc(size_t size) {
    double* ptr;
    if (posix_memalign(reinterpret_cast<void**>(&ptr), 64, size * sizeof(double)) != 0) {
        return std::unique_ptr<double[], decltype(&std::free)>(nullptr, std::free);
    }
    return std::unique_ptr<double[], decltype(&std::free)>(ptr, std::free);
}

// Funkcja obliczająca odległość euklidesową bez użycia AVX-512
double euclidean_distance_no_avx(const double* a, const double* b, size_t length) {
    double sum = 0.0;
    for (size_t i = 0; i < length; ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }
    return std::sqrt(sum);
}

// [[Rcpp::export]]
NumericMatrix dist_no_avx_parallel(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    auto aligned_Z = aligned_alloc(n * dim);

    if (!aligned_Z) {
        stop("Allocation failed");
    }

    // Kopiowanie danych do wyrównanych buforów
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            aligned_Z[i * dim + j] = Z(i, j);
        }
    }
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double* vec_i = aligned_Z.get() + i * dim;
            const double* vec_j = aligned_Z.get() + j * dim;

            double dist = euclidean_distance_no_avx(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;
    }

    return DZ;
}

// [[Rcpp::export]]
NumericMatrix dist_no_avx_no_parallel(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    auto aligned_Z = aligned_alloc(n * dim);

    if (!aligned_Z) {
        stop("Allocation failed");
    }

    // Kopiowanie danych do wyrównanych buforów
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            aligned_Z[i * dim + j] = Z(i, j);
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double* vec_i = aligned_Z.get() + i * dim;
            const double* vec_j = aligned_Z.get() + j * dim;

            double dist = euclidean_distance_no_avx(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;
    }

    return DZ;
}