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
**TODO**

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
- **Suffix array querying** - implementation in *O*(*n* + *log m*) adapted from prof. Szczurek's notes [\[LINK\]](https://www.mimuw.edu.pl/~szczurek/TSG2/04_suffix_arrays.pdf)
- **All other code** comes from classes or is a product of own work
