# Assignment 1: read mapping

## Solution
Given the reference string initialise `Index` by computing its suffix array using provided Karkkainen-Sanders implementation (`src/KarkSand.py`). Then in a loop, for each read:
1. Compute a set of 0- and 1-edit-distance shingles from the read string. For each shingle memorize its sequence, offset from the beginning of read sequence and range of estimated occurence positions in reference.
    - Length and amount of shingles depends on different parameters, such as length of read, expected error rate, or step in indexes between neighbouring ones,
    - First compute exact shingles from the read,
    - Then for each shingle compute its 1-neighbourhood.
2. For each shingle in the set:
    - Query suffix array for its position,
    - Compute its range of estimated occurence positions in reference,
    - Add `+1` and `-1` to a cumulative occurence map at range's start and stop positions accordingly.
3. For each entry of `(index, count)` in the sorted cumulative occurence map:
    - Add `count` to cumulative count,
    - If cumulative count is the greatest so far remember the `index` as most likely matching occurence.
4. With the index of most likely read matching occurence in reference compute the `KEditDP` class to find best edit distance and its starting and ending positions. Implementation was optimised to compute the matrix diagonally, allowing to take advantage of numpy vectorisation.

## Complexity
| Variable | Definition                                      |
| -------- | ----------------------------------------------- |
| $P$      | Length of searched pattern (read)               |
| $R$      | Length of reference string                      |
| $S$      | Length of one shingle                           |
| $N$      | Number of shingles generated per pattern (read) |
| $D$      | Distance between shingles' starting positions   |
| $E$      | Chosen error rate                               |

**Note:** $D$ is set to be around $10$ in the presented implementation to match the one second per read requirement.

The following estimates hold for the presented implementation:
- $S \approx 2 \div E$
- $N \approx P \cdot S \cdot \frac{6}{D}$
    - $6$ comes from generating $1$-neighbourhood. For one shingle of length $S$ there are $6$ possible operations at each position:
        - $1$ substitution
        - $1$ deletion
        - $4$ insertions

The key time complexity factors are:
- $O(R)$ - Initialising the index. However, it only has to be done once for the given reference string and can serve any number of reads.
- $O(S + log(R))$ - Querying the suffix array for occurences of one shingle. (Note that this is amorthised complexity. For it to be a fixed complexity we would have to sacrifice memory for storing LCP arrays, but we cannot do that with existing  $1\,GB$ memory limit).
- $O(N \cdot log (N))$ - Sorting the occurence range map.
- $O(\cdot P^2)$ - Computing the `KEditDP`.

This allows us to finally write down the complexity for processing a single pattern (read), assuming the index was already initialised:
$$ O(N \cdot (S + log(R)) + N \cdot log(N) + P^2) $$
what can be simplified to:
$$ O(N \cdot (S + log(R \cdot N)) + P^2) $$
If we were to substitute the values from the implementation $S \approx 20$, $N \approx 1.3 \cdot 10^4$ while assuming the maximum reference length $R = 2 \cdot 10^7$ we can see that:
- $O(N \cdot (S + log(R \cdot N))) \approx 8 \cdot 10^6$
- $O(P^2) \approx 10^6$

What means that we have implemented the solution choosing the right values of $S$ and $N$ to achieve high accuracy without exceeding the unavoidable computational overhead of $O(P^2)$

## Plots
The script provides ability to plot the result of shingle occurence mapping over the reference string. Seeing the plots is a great way of verifying the script's accuracy and understanding its work principles. To generate plots uncomment the line
```
# self.__plot_occurence_map(occurence_map, read_id)
```
inside `query()` method of `Index` class. Plots will be saved to `tests/plots` directory with names matching `read_id`'s.

## Tests
Tests are structured in two subdirectories, `sample` and `large`. They were included together with the assignment. `sample` tests are good for a quick, possibly inaccurate verification of the program. `large` tests are great for experimenting with different parameter values, as well as testing the script against some more difficult examples.

Two scripts are provided for convenient testing of the script.

1. `run.sh` is used for running the script on a selected test directory. The results will be saved to `output` directory.
```
./run.sh <test_directory_path>
```

2. `statistics.sh` is used for comparing the `run.sh` output with the expected output, verifying it against the specified thresholds. It counts read accuracy, number of inaccuracies, as well as runtime, and computes the estimated number of points awarded for the current version of the script. Also `metrics.txt` file will be generated and saved in the test directory for future inspection.
```
./statistics.sh <test_directory_path>
```

## Bibliography
- **Karkkainen-Sanders algorithm** - implementation adapted from Google archive [\[LINK\]](https://code.google.com/archive/p/pysuffix/)
- **Suffix array querying** - implementation in $O(n + log(m))$ adapted from prof. Szczurek's notes [\[LINK\]](https://www.mimuw.edu.pl/~szczurek/TSG2/04_suffix_arrays.pdf)
- **All other code** comes from classes or is a product of own work
