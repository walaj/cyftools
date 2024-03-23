# Cyftools
A command-line interface for wrangling data from CyCIF images

## Build

### Clone
Clone the project recursively (it has submodules):

```
git clone --recursive git@github.com:walaj/cyftools
```

### Install Dependencies
Install the following dependencies:

- libjpeg
- hdf5
- eigen3
- libomp
- cgal
- boost
- mlpack
- libbdplus
- armadillo

On OSX this can be done with:

```
brew install libjpeg hdf5 eigen libomp cgal boost mlpack libbdplus armadillo
```

#### supervised-lda library
You will also need `supervised-lda`, which needs to be in `$HOME/git/supervised-lda`, so run:

```
mkdir -p $HOME/git
git clone git@github.com:angeloskath/supervised-lda $HOME/git/supervised-lda
```

Then follow the instructions [here](https://github.com/angeloskath/supervised-lda/blob/master/CMakeLists.txt) to build it from source (might need to install cmake first as well).

On MacBook Arm architectures the `make check` on supervised-lda might fail with stuff like:

```
test_approximate_supervised_expectation_step.cpp:166:5: error: reference to 'VectorX' is ambiguous`
fatal error: too many errors emitted, stopping now [-ferror-limit=]
```

Running `sudo make install` and then `cp build/libldaplusplus.dylib lib/` got it to work anyway - `cyftools` expects the library to be in `$HOME/git/supervised-lda/lib/libldaplusplus.dylib` rather than `build`.

### Build from Source
From root of repo:
```
cd src
make all
```
