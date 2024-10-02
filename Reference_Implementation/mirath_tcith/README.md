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
cmake -DCMAKE_BUILD_TYPE=<BUILD_TYPE> -DOPT_VAR=<VAR_TYPE> -DOPT_SEC=<SEC_TYPE> -DOPT_LEN=<LEN_TYPE> -B build
```

where `<BUILD_TYPE>` is either **Release** or **Debug**, `<VAR_TYPE>` is either **a** or **b**, `<SEC_TYPE>` is either **I**, **III**, and **V**, and `<LEN_TYPE>` is either **SHORT** or **FAST**.

Optionally, add personal preference to CMAKE such as `-DCMAKE_MAKE_PROGRAM=path/to/ninja -G Ninja`

## Example

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -B build
cmake -DCMAKE_BUILD_TYPE=Debug -DOPT_VAR=a -DOPT_SEC=I -DOPT_LEN=FAST -B build
cmake -DCMAKE_BUILD_TYPE=Debug -DOPT_VAR=a -DOPT_SEC=III -DOPT_LEN=FAST -B build
cmake -DCMAKE_BUILD_TYPE=Debug -DOPT_VAR=a -DOPT_SEC=V -DOPT_LEN=FAST -B build
```

Jump into buiild folder and build all target libraries for all supported parameter sets:

```bash
cd build
make -j8
```