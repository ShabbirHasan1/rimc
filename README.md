# rsimc

A rust simulator for Monte Carlo systems. Initially built of off the Ising model, this is both an executable and a library that can be used in other systems.

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

### Configure simulation

*TBD*

Ideally, the input is intended to use a [toml](https://toml.io/en/) to set up the parameters of the simulation. An example of this would be,

```
cargo run - example/config.toml
```

```
cd target/release
./rsimc $PRJ_ROOT/example/config.toml
```


