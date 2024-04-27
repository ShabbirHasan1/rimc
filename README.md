# rs

A rust comparison for Ising model Monte Carlo systems. I am comparing: rust, python and C++ work to develop a system. The key points here are,

* Code readability
* Testing
* Ecosystem
* Confidence in the code

I found that while rust in new to me, building the system was quite pleasant. Due to the way rust has all the built-in features for testing and compile time checks I was surprised (in shock?) that the code ran the first time. Python had several run time bugs that needed to be fixed when transcribing the code. C++ was AI generated based on the completed python code, which was an interesting choice to see how I could quickly gain speed from a python implementation. While this code base is quite small and easy to review, I would not trust the C++ AI transcribe of python code for larger code bases that require extensive review.

Overall, this proved to be quite a useful exercise in development in code and I believe that all the above points were easily reached with rust.

## Build

```
cargo build
```

### Run

You can either use cargo to run or run the executable natively

```
cargo run -- -d 20
```

or 

```
./target/build/debug/rsimc -d 20
```

## Benchmark

We compare a slightly optimized python implementation of the same algorithm to both a debug and release version of the rust equivalent. You can see below that we get about 3x or 54x the speed improvement when using *debug* or *release* mode respectively.

![Benchmark](./assets/benchmark.png)


A *lazy* C++ implementation was added, which *should* be comparable in speed to the rust equivalent. We used ChatGPT4 to rewrite the python code implementation to the `-std=c++11 -O3` and keep external packages to a minimum. Note that this does not take into account any development hurdles. As the rust implementation was only compiled once, with only one *quasi bug* due to how the modulus `%` operator was assumed to work, which really was a remainder operator. With rust, we needed `a.rem_euclid(b)` for the equivalent `a % b` operation.

Rust seems to perform equally well or up to 2.5x faster when compared to C++. The peak difference occurs at a system size of N=1000.
