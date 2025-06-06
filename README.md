# distinct-kmers

This is a small (but fairly optimized) tool that counts the number of distinct k-mers in a sequence.
It supports (possibly compressed and multi-line) fasta files and k <= 32.

To use it, simply clone this repository and run
```sh
cargo r -r -- -k <K> -i <FASTA/FASTQ/@FOF> --fastq
```

In case of File of File, on line per input file in the fof file. Inputs must all be either fastq of fasta. 

Alternatively, you can install the current version locally with
```sh
RUSTFLAGS="-C target-cpu=native" cargo install -f --git https://github.com/pierrepeterlongo/distinct-kmers.git
```
This will create a binary named `distinct-kmers` and add it to your `PATH`.

## Optimizations

- [x] native build with [LTO](https://nnethercote.github.io/perf-book/build-configuration.html#link-time-optimization)
- [x] compute local counts using super-k-mer buckets
- [x] SIMD computation of super-k-mers using [`simd-minimizers`](https://crates.io/crates/simd-minimizers)
- [x] transparent decompression using [`niffler`](https://crates.io/crates/niffler)
- [x] parallel FASTA parsing using [`seq_io_parallel`](https://crates.io/crates/seq_io_parallel)
- [x] parallel local counting using [`rayon`](https://crates.io/crates/rayon)
- [x] faster hash using [`FxHash`](https://crates.io/crates/rustc-hash)

Failed attempts:
- [ ] using [`boxcar`](https://crates.io/crates/boxcar) for concurrent buckets (vectors with mutex ended up being faster)
- [ ] using [`mimalloc`](https://crates.io/crates/mimalloc) (no improvement over the default allocator)

## Acknowledgement

This tool was inspired by https://github.com/pierrepeterlongo/unique_kmer_counter.
My goal here was to create a simple yet optimized implementation for the same task.
