# Poprawa wydajności wielowymiarowego testu statystycznego typu Wilcoxona–Manna–Whitneya
### Błażej Pałkus 
### Wydział Matematyki, Informatyki i Mechaniki Uniwersytetu Warszawskiego

Repozytorium zawiera kod źródłowy oraz pliki danych związane z pracą magisterską "Poprawa wydajności wielowymiarowego testu statystycznego typu Wilcoxona–Manna–Whitneya". Pełny tekst pracy można znaleźć [tutaj](https://apd.uw.edu.pl/diplomas/231089/).

## Struktura repozytorium

- `setup.R`: Skrypt instalujący wymagane pakiety R oraz źródła C++.
- `rank_distance_original.R`: Kod oryginalnej funkcji autorów pierwotnej pracy.
- `rank_distance_upgraded.R`: Kod funkcji zoptymalizowanej bez akceleracji sprzętowej.
- `rank_distance_upgraded_cpp.R`: Kod funkcji zoptymalizowanej z akceleracją sprzętową realizowaną za pomocą kodu C++.
- `distance_avx512.cpp`: Kod w C++ obliczający odległości euklidesowe z użyciem AVX-512.
- `distance_avx512_float.cpp`: Kod w C++ obliczający odległości euklidesowe z użyciem AVX-512 dla typów float.
- `distance_no_avx512.cpp`: Kod w C++ obliczający odległości euklidesowe bez użycia AVX-512.
- `distance_no_avx512_float.cpp`: Kod w C++ obliczający odległości euklidesowych bez użycia AVX-512 dla typów float.
- `avx_support.cpp`: Kod w C++ sprawdzający wsparcie dla AVX-512.
- Pliki danych eksperymentalnych w formacie `.csv`.

## Wymagania

- R (wersja >= 3.6.0)
- Kompilator C++ obsługujący standard C++11 lub wyższy
- Biblioteki C++: OpenMP (oraz ew. AVX-512)

## Instrukcja instalacji

Aby uruchomić kod, należy wykonać poniższe kroki:

1. Sklonować repozytorium:
    ```bash
    git clone https://github.com/BlazejUW/multiWMW
    cd multiWMW
    ```

2. Uruchomić skrypt `setup.R`, który zainstaluje wymagane pakiety i skompiluje oraz wczyta pliki C++:
    ```r
    source("setup.R")
    ```

## Opis funkcji

### rank_distance_original

Oryginalna funkcja autorów pracy "A generalized Wilcoxon–Mann–Whitney type test for multivariate data through pairwise distance". U autorów pracy występuję pod nazwą `Rank_distance2`, w mojej pracy występuję pod zmienioną nazwą `rank_distance_original`. Znajduje się w pliku `rank_distance_original.R`. Oblicza statystykę na podstawie odległości euklidesowych pomiędzy punktami w przestrzeni wielowymiarowej.

### rank_distance_upgraded

Zoptymalizowana wersja funkcji bez akceleracji sprzętowej, poprawia wydajność obliczania statystyki testowej poprzez zmniejszenie złożoności obliczeniowej z `O(k^3)` do `O(k^2 log(k))`. Znajduje się w pliku `rank_distance_upgraded.R`.

### rank_distance_upgraded_cpp

Zoptymalizowana wersja funkcji obliczającej statystykę testową z akceleracją sprzętową realizowaną za pomocą kodu C++, wykorzystująca wsparcie dla AVX-512 oraz możliwość równoległego przetwarzania. Znajduje się w pliku `rank_distance_upgraded_cpp.R`.

#### Użycie funkcji rank_distance_upgraded_cpp

Funkcja `rank_distance_upgraded_cpp` przyjmuje dodatkowe flagi jako argumenty:
- `use_avx512`: Flaga wskazująca, czy używać instrukcji AVX-512 (domyślnie TRUE). Jeśli procesor nie obsługuje AVX-512, program to sprawdzi i automatycznie zmieni na wersję bez użycia AVX-512.
- `use_parallel`: Flaga wskazująca, czy używać przetwarzania równoległego (domyślnie TRUE).
- `use_float`: Flaga wskazująca, czy używać typów float zamiast double (domyślnie FALSE).

## Pliki C++

### distance_avx512.cpp

Plik `distance_avx512.cpp` zawiera funkcje do obliczania macierzy odległości euklidesowych z użyciem instrukcji AVX-512 dla typów double.

### distance_avx512_float.cpp

Plik `distance_avx512_float.cpp` zawiera funkcje do obliczania macierzy odległości euklidesowych z użyciem instrukcji AVX-512 dla typów float.

### distance_no_avx512.cpp

Plik `distance_no_avx512.cpp` zawiera funkcje do obliczania macierzy odległości euklidesowych bez użycia instrukcji AVX-512, przeznaczony dla maszyn, które nie obsługują AVX-512.

### distance_no_avx512_float.cpp

Plik `distance_no_avx512_float.cpp` zawiera funkcje do obliczania macierzy odległości euklidesowych dla typów float bez użycia instrukcji AVX-512.

### avx_support.cpp

Plik `avx_support.cpp` zawiera funkcję sprawdzającą wsparcie dla instrukcji AVX-512 na danej maszynie.

## Pliki z wynikami eksperymentów

Pliki z wynikami eksperymentów znajdują się w repozytorium i kończą się rozszerzeniem `.csv`. Szczegółowy opis eksperymentów i wyników znajduje się w pracy magisterskiej.

## Kontakt

W przypadku pytań dotyczących mojej pracy magisterskiej, proszę o kontakt poprzez otwarcie issue w repozytorium GitHub lub bezpośrednio poprzez email (blazej.palkus@gmail.com).
