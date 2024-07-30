#include <omp.h> // OpenMP header
#include <Rcpp.h>
#include <immintrin.h> // AVX-512 header
#include <cmath>
#include <memory>
#include <iostream>

using namespace Rcpp;

// Funkcja alokująca wyrównaną pamięć
std::unique_ptr<float[], decltype(&_mm_free)> aligned_alloc_float(size_t size) {
    return std::unique_ptr<float[], decltype(&_mm_free)>(
        (float*)_mm_malloc(size * sizeof(float), 64), _mm_free);
}

// Funkcja obliczająca odległość euklidesową z użyciem AVX-512 dla float
float euclidean_distance_float(const float* a, const float* b, size_t length) {
    __m512 vec_sum = _mm512_setzero_ps();
    size_t i = 0;

    for (; i + 16 <= length; i += 16) {
        __m512 vec_a = _mm512_loadu_ps(a + i);
        __m512 vec_b = _mm512_loadu_ps(b + i);
        __m512 diff = _mm512_sub_ps(vec_a, vec_b);
        __m512 sq_diff = _mm512_mul_ps(diff, diff);
        vec_sum = _mm512_add_ps(vec_sum, sq_diff);
    }

    alignas(64) float buffer[16];
    _mm512_storeu_ps(buffer, vec_sum);
    float sum = 0.0f;
    for (int j = 0; j < 16; ++j) {
        sum += buffer[j];
    }

    for (; i < length; ++i) {
        float diff = a[i] - b[i];
        sum += diff * diff;
    }

    __m512 vec_sum_sqrt = _mm512_sqrt_ps(_mm512_set1_ps(sum));
    _mm512_storeu_ps(buffer, vec_sum_sqrt);
    return buffer[0];
}

// Funkcja obliczająca macierz odległości euklidesowych z użyciem float
// [[Rcpp::export]]
NumericMatrix dist_avx512_float(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    // Alokacja wyrównanych buforów na całą macierz
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

            float dist = euclidean_distance_float(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;  // Dystans do samego siebie jest zawsze 0
    }

    return DZ;
}

// Funkcja obliczająca macierz odległości euklidesowych z użyciem float
// [[Rcpp::export]]
NumericMatrix dist_avx512_parallel_float(const NumericMatrix& Z) {
    size_t n = Z.nrow();
    size_t dim = Z.ncol();
    NumericMatrix DZ(n, n);

    // Alokacja wyrównanych buforów na całą macierz
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

            float dist = euclidean_distance_float(vec_i, vec_j, dim);
            DZ(i, j) = dist;
            DZ(j, i) = dist;
        }
        DZ(i, i) = 0.0;  // Dystans do samego siebie jest zawsze 0
    }

    return DZ;
}