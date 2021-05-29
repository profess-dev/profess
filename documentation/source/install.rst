*******
Install
*******

The simplest approach is to first install the various dependencies into a conda environment. ::

    conda create --name profess cmake python=3 numpy scipy
    conda activate profess
    conda install -c conda-forge libxc armadillo fftw

Next, clone the profess repo, using ``--recursive`` to obtain necessary submodules. ::

    git clone --recursive https://github.com/profess-dev/profess
    cd profess

Then, build ``profess`` with cmake, replacing <conda-dir> with the path to your conda root directory. ::

    mkdir build; cd build
    cmake -DLIBXC_INCLUDE_DIR=<conda-dir>/envs/profess/include -DLIBXC_LIBRARY_DIR=<conda-dir>/envs/profess/lib ..
    make
    export PYTHONPATH=$PYTHONPATH:$PWD

Finally, run the tests. ::

    cd ../test
    python -m unittest
