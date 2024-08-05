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
- `reducing_z_ratio_without_redistribution.py`: Skrypt Python do eksperymentów nad zmniejszeniem liczności próbki referencyjnej \(Z\) bez redystrybucji.
- `reducing_z_ratio_with_redistribution.py`: Skrypt Python do eksperymentów nad zmniejszeniem liczności próbki referencyjnej \(Z\) z redystrybucją do próbek \(X\) i \(Y\).
- Pliki zawierające wyniki eksperymentów w formacie `.csv`.

## Wymagania

- R (wersja >= 3.6.0)
- Kompilator C++ obsługujący standard C++11 lub wyższy
- Biblioteki C++: OpenMP (oraz ew. AVX-512)
- Python (wersja >= 3.6)
- Pakiety Python: numpy, scipy, pandas, scikit-learn, joblib, multiprocessing

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

3. Zainstalować wymagane pakiety Python:
    ```bash
    pip install numpy scipy pandas scikit-learn joblib
    ```

## Opis funkcji

### rank_distance_original

Oryginalna funkcja autorów pracy "A generalized Wilcoxon–Mann–Whitney type test for multivariate data through pairwise distance". U autorów pracy występuje pod nazwą `Rank_distance2`, w mojej pracy występuje pod zmienioną nazwą `rank_distance_original`. Znajduje się w pliku `rank_distance_original.R`. Oblicza statystykę na podstawie odległości euklidesowych pomiędzy punktami w przestrzeni wielowymiarowej.

### rank_distance_upgraded

Zoptymalizowana wersja funkcji bez akceleracji sprzętowej, poprawia wydajność obliczania statystyki testowej poprzez zmniejszenie złożoności obliczeniowej z `O(k^3)` do `O(k^2 log(k))`. Znajduje się w pliku `rank_distance_upgraded.R`.

### rank_distance_upgraded_cpp

Zoptymalizowana wersja funkcji obliczającej statystykę testową z akceleracją sprzętową realizowaną za pomocą kodu C++, wykorzystująca wsparcie dla AVX-512 oraz możliwość równoległego przetwarzania. Znajduje się w pliku `rank_distance_upgraded_cpp.R`.

#### Użycie funkcji rank_distance_upgraded_cpp

Funkcja `rank_distance_upgraded_cpp` przyjmuje dodatkowe flagi jako argumenty:
- `use_avx512`: Flaga wskazująca, czy używać instrukcji AVX-512 (domyślnie TRUE). Jeśli procesor nie obsługuje AVX-512, program to sprawdzi i automatycznie zmieni na wersję bez użycia AVX-512.
- `use_parallel`: Flaga wskazująca, czy używać przetwarzania równoległego (domyślnie TRUE).
- `use_float`: Flaga wskazująca, czy używać typów float zamiast double (domyślnie FALSE).

## Użycie funkcji w R
Aby uruchomić każdą z funkcji, można użyć przykładowych danych:
```r
source("rank_distance_original.R") # Wczytanie funkcji rank_distance_original
source("rank_distance_upgraded.R") # Wczytanie funkcji rank_distance_upgraded
source("rank_distance_upgraded_cpp.R") # Wczytanie funkcji rank_distance_upgraded_cpp
X <- matrix(rnorm(100), ncol=5)
Y <- matrix(rnorm(100), ncol=5)
Z <- matrix(rnorm(200), ncol=5)
result_original <- rank_distance_original(X, Y, Z)
result_upgraded <- rank_distance_upgraded(X, Y, Z)
result_upgraded_cpp <- rank_distance_upgraded_cpp(X, Y, Z)
```
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

## Skrypty Python

### reducing_z_ratio_without_redistribution.py

Skrypt `reducing_z_ratio_without_redistribution.py` służy do eksperymentów nad zmniejszeniem liczności próbki referencyjnej \(Z\) bez redystrybucji. Skrypt generuje próbki \(X\), \(Y\) i \(Z\), a następnie zmniejsza liczność \(Z\) do różnych poziomów (75%, 60%, 50%, 37%, 25%, 10%, 5%, 1%) i analizuje wpływ tych zmian na p-wartości. Skrypt wykorzystuje metodę bootstrapu, która jest uruchamiana 500 razy, aby oszacować rozkład statystyki testowej i obliczyć p-wartości.

### reducing_z_ratio_with_redistribution.py

Skrypt `reducing_z_ratio_with_redistribution.py` służy do eksperymentów nad zmniejszeniem liczności próbki referencyjnej \(Z\) z redystrybucją do próbek \(X\) i \(Y\). Skrypt generuje próbki \(X\), \(Y\) i \(Z\), następnie zmniejsza liczność \(Z\) do różnych poziomów, redystrybuuje elementy \(Z\) do \(X\) i \(Y\) oraz analizuje wpływ tych zmian na p-wartości i czas działania algorytmu. Skrypt również wykorzystuje metodę bootstrapu, która jest uruchamiana 500 razy, aby oszacować rozkład statystyki testowej i obliczyć p-wartości.

## Uruchamianie skryptów Python
Aby uruchomić skrypty Python, należy uruchomić poniższe komendy:
```bash
python reducing_z_ratio_without_redistribution.py
python reducing_z_ratio_with_redistribution.py
```
Dane do eksperymentów są generowane losowo i automatycznie. Zapisują się w odpowiednich plikach.
## Pliki z wynikami eksperymentów

Pliki z wynikami eksperymentów znajdują się w repozytorium i kończą się rozszerzeniem `.csv`. Szczegółowy opis eksperymentów i wyników znajduje się w pracy magisterskiej.

## Kontakt

W przypadku pytań dotyczących mojej pracy magisterskiej, proszę o kontakt poprzez otwarcie issue w repozytorium GitHub lub bezpośrednio poprzez email (blazej.palkus@gmail.com).
