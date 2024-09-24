# Mirath

## Description

...

### Current Supported Parameter Sets

We support the following parameter sets (Ia, IIIa, and Va).

## Compilation

Compilation for Linux and MacOS using CMAKE

### cmake configuration

First configure the build directory with

```bash
cmake -DCMAKE_BUILD_TYPE=<BUILD_TYPE> -DOPT_SET=<OPT_TYPE> -B build
```

where `<BUILD_TYPE>` is either **Release** or **Debug**.
Additionally, `<OPT_TYPE>` is either **Ia**, **IIIa** or **Va**; by default it assumes **Ia**.

Optionally, add personal preference to CMAKE such as `-DCMAKE_MAKE_PROGRAM=path/to/ninja -G Ninja`

## Example

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -B build
cmake -DCMAKE_BUILD_TYPE=Debug -DP_SET=Ia -B build
cmake -DCMAKE_BUILD_TYPE=Debug -DP_SET=IIIa -B build
cmake -DCMAKE_BUILD_TYPE=Debug -DP_SET=Va -B build
```

Jump into buiild folder and build all target libraries for all supported primes:

```bash
cd build
make -j8
```