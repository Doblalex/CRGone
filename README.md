# Pace Challenge 2024 Solver CRGone

Requirements: [gmp](https://gmplib.org/#DOWNLOAD)

```
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make -j4
./olcmexact < input.gr
```

Do not use too many threads for `make` as scip will get stuck.