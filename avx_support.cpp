#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// Funkcja sprawdzająca wsparcie dla AVX-512
// [[Rcpp::export]]
bool check_avx512_support() {
#ifdef __AVX512F__
    return true;
#else
    return false;
#endif
}