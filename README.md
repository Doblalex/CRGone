```
git submodule update --init --recursive
mkdir build
cd buid
cmake ..
make -j4
./olcmexact < input.gr
```

Do use too many threads for `make` as scip will get stuck.