# Cyftools
A command-line interface for wrangling CyCIF images

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

On OSX this can be done with:

```
brew install libjpeg hdf5 eigen libomp cgal boost mlpack
```

You will also need `supervised-lda`, which needs to be in `$HOME/git/supervised-lda`, so run:

```
mkdir -p $HOME/git
git clone git@github.com:angeloskath/supervised-lda $HOME/git/supervised-lda
```

Then follow the instructions [here](https://github.com/angeloskath/supervised-lda/blob/master/CMakeLists.txt) to build it from source (might need to install cmake first as well). 

### Build from Source
From root of repo:
```
cd src
make all
```
