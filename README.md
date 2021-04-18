# Approximate Shannon Entropy

A rust library to calculate the approximate Pearson Median Skew of a slice.

The median algorithm used is a O(1) size but O(nk) time complexity with worst case O(n^2), due to the immutable slice input restriction.

Usable on `no_std` due to use of approximate square root from [micromath](https://github.com/tarcieri/micromath).

# Usage

Add this to your Cargo.toml
```
[dependencies]
approx_pearson_skew = "0.1.0"
```

# Examples

```
$ cargo build --example stdin_skew
$ echo 1234 | ./target/debug/examples/stdin_skew
```