#include <omp.h> // OpenMP header
#include <Rcpp.h>
#include <immintrin.h> // AVX-512 header
#include <cmath>
#include <memory>
#include <iostream>
using namespace Rcpp;

// Funkcja alokująca wyrównaną pamięć
std::unique_ptr<double[], decltype(&_mm_free)> aligned_alloc(size_t size) {
    return std::unique_ptr<double[], decltype(&_mm_free)>(
        (double*)_mm_malloc(size * sizeof(double), 64), _mm_free);
}

// Funkcja obliczająca odległość euklidesową z użyciem AVX-512
double euclidean_distance(const double* a, const double* b, size_t length) {
    __m512d vec_sum = _mm512_setzero_pd();
    size_t i = 0;

    for (; i + 8 <= length; i += 8) {
        __m512d vec_a = _mm512_loadu_pd(a + i);
        __m512d vec_b = _mm512_loadu_pd(b + i);
        __m512d diff = _mm512_sub_pd(vec_a, vec_b);
        __m512d sq_diff = _mm512_mul_pd(diff, diff);
        vec_sum = _mm512_add_pd(vec_sum, sq_diff);
    }

    alignas(64) double buffer[8];
    _mm512_storeu_pd(buffer, vec_sum);
    double sum = 0.0;
    for (int j = 0; j < 8; ++j) {
        sum += buffer[j];
    }

    for (; i < length; ++i) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }

    __m512d vec_sum_sqrt = _mm512_sqrt_pd(_mm512_set1_pd(sum));
    _mm512_storeu_pd(buffer, vec_sum_sqrt);
    return buffer[0];
}

// Funkcja obliczająca macierz odległości euklidesowych
// [[Rcpp::export]]
NumericMatrix dist_avx512(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    // Alokacja wyrównanych buforów na całą macierz
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

            double dist = euclidean_distance(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;  // Dystans do samego siebie jest zawsze 0
    }

    return DZ;
}

// [[Rcpp::export]]
NumericMatrix dist_avx512_parallel(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    // Alokacja wyrównanych buforów na całą macierz
    auto aligned_Z = aligned_alloc(n * dim);

    if (!aligned_Z) {
        stop("Allocation failed");
    }

    // Kopiowanie danych do wyrównanych buforów
    #pragma omp parallel for
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

            double dist = euclidean_distance(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;  // Dystans do samego siebie jest zawsze 0
    }

    return DZ;
}