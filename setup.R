# Ustawienie globalnego repozytorium CRAN
options(repos = c(CRAN = "https://cran.r-project.org"))

install_and_load <- function(packages) {
  for (package in packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package, dependencies = TRUE)
      library(package, character.only = TRUE)
    }
  }
}

# Przykład użycia
install_and_load(c("testthat", "ggplot2", "microbenchmark", "Rcpp", "dplyr", "MASS", "mnormt", "snow", "mvtnorm"))

sourceCpp("avx_support.cpp")
sourceCpp("distance_no_avx512.cpp")
sourceCpp("distance_no_avx512_float.cpp")
if(check_avx512_support()){
    Sys.setenv("PKG_CXXFLAGS"="-O2 -mavx512f -mavx512bw -mavx512cd -mavx512dq -mavx512vl")
    sourceCpp("distance_avx512.cpp")
    sourceCpp("distance_avx512_float.cpp")
}
