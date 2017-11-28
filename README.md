## part_conc_spt_abc

Nanoparticle concentration measurements using single particle tracking and
Approximate Bayesian Computation.


### Instructions: Running the CPU version

The code is developed using Julia 0.6.0. To run the analysis, just run src/cpu/run_part_conc_spt_abc.jl in the following way

julia -p <number_of_cores> <program> <input>

where the parameters (all mandatory) are

- `<number_of_cores>` : number of CPU cores (threads) to utilize
- `<program>` : path to run_part_conc_spt_abc.jl
- `<input>` : input file

The file test/test_part_conc_spt_abc.jl generates a simulated example data set
and calls the analysis program.

### Instructions: Building \& Running the GPU version

**Building the GPU version:**

First generate the makefiles using [premake](https://premake.github.io/):
```
$ external/premake5 gmake
```
Note that only make-based builds have been tested. The project definition in
`premake5.lua` includes raw `nvcc` command lines which may or may not need to
be adjusted for other build types ("actions") supported by premake.

Next, run make:
```
$ make -j4
```
This builds the debug version. An optimized release version can be built via
```
$ make -j4 config=release_native
```

This should produce the main binary (`accel-<config>.exe`) in the `bin/`
directory.

**Running the GPU version:**

Run `bin/accel-<config>.exe --help` for an overview of command line options.
The most important options are:

 - `-i <input>` : select input file
 - `-o <output>` : select output file (will be overwritten!)
 - `-g <gpuspec>` : select GPU configuration (see below)

Other command line options can typically be omitted (which will use their
default values).

If no output file is specified via `-o`, the default output file path from the
input will be used.

The `<gpuspec>` is used to select GPUs. It is a comma-separated list (no
spaces!) of `<GPU number>/<number of queues>`. The default is `-g 1/30`, which
uses the first available GPU with 30 asynchronous queues. For a system with two
identical GPUs one might use `-g 1/30,2/30`. Work is distributed evently across
the queues, so simple manual "load-balancing" can be performed by lowering the
number of queues (e.g., for a system with a GTX1060 + a GTX960, `-g 1/30,2/21`
turned out to be a good choice).

Available CUDA GPUs can be listed by using the `--list-devices` command line
option.
