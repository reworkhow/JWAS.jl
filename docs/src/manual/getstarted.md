# Get Started

## Installation

To install julia, please go to the [offical Julia website](https://julialang.org/downloads/).
Please see [platform specific instructions](https://julialang.org/downloads/platform.html)
if you have trouble installing Julia.

To install the package, use the following command inside the Julia REPL (or IJulia Notebook):
```julia
Pkg.add("JWAS")
```

To load the JWAS package, use the following command inside the Julia REPL (or IJulia Notebook):

```julia
using JWAS
```

The command `Pkg.add("JWAS")` will add the registered official JWAS package and dependencies.

To use the latest/beta features under development, run `Pkg.checkout("JWAS")` to get the
newest unofficial JWAS. Run `Pkg.free("JWAS")` to go back to the offical one.

###### IJulia notebook

If you prefer “reproducible research”, an interactive Jupyter notebook interface is available
for Julia (and therefore JWAS). The Jupyter notebook is an open-source web application for creating
and sharing documents that contain live code, equations, visualizations and explanatory text.
To install IJulia, please go to [IJulia](https://github.com/JuliaLang/IJulia.jl).

###### Jupyter-IJulia notebooks via Docker

Docker provides a straightforward way to install Jupyter-IJulia notebooks with JWAS.

- Install Docker from [here](https://docs.docker.com/install/) for your platform.

- From a terminal (on Mac or Linux), run the command:

```bash
docker run -it --rm -p 8888:8888 qtlrocks/jwas-docker
```

This will start a Jupyter-IJulia Notebook server listening for HTTP connections on port 8888 with a randomly generated authentication
token. Documentation for JWAS can be accessed from the notebook: "JWAS_notebooks/index.ipynb".

The directories and files created within the Docker container will be lost when the container is stopped. To save your work
on the host machine, a directory on the host machine can be mounted as a folder in the container with the command:

```bash
docker run -it --rm -p 8888:8888 -v path_to_folder_on_host:/home/jovyan/folder_in_container qtlrocks/jwas-docker
```

where `path_to_folder_on_host` is the path to the folder that you want to have access to from within the container, and  
`folder_in_container` is the name of the folder in the container. For example, the Docker command

```bash
docker run -it --rm -p 8888:8888 -v /Users/rohan:/home/jovyan/rohan qtlrocks/jwas-docker
```

creates a Docker container with the folder `rohan` with the contents of `/Users/rohan` of the host machine. Files and
directories that are in the folder `/home/jovyan/rohan` will not be lost when the container is stopped.  



###### Standalone application

!!! note "standalone application (no installation required)"

    A fully self-contained application for JWAS (no installation required) will come out this year.


## Access documentation

To show the basic information (README file) of JWAS in REPL or IJulia notebook using `?JWAS`
and press enter.

For help on a specific function, type ? followed by its name, e.g. `?runMCMC` and press enter
in REPL or IJulia notebook.

!!! warning

    Please load the JWAS package at first.


The full documentation is available [here](http://reworkhow.github.io/JWAS.jl/latest/index.html).

## Run your analysis

There are several ways to run your analysis.

(1) The easiest way to run analysis in Julia is by starting **an interactive session (REPL)** by double-clicking the Julia
executable or running julia from the command line (e.g., terminal) as

```julia
julia> 1+2
3

julia> 3*4
12
```

To evaluate code written in a file `script.jl` in REPL, write and run

```julia
julia> include("script.jl").
```
To exit the interactive session, type `^D` – the control key together with the d key or type `quit()`.

(2) To run code in a file non-interactively from **the command line** (e.g.,termial), you can give it as the first argument to the `julia` command:

```bash
julia script.jl
```

If you want to pass arguments to your script, run it as
```bash
julia script.jl arg1 arg2
```
where arguments `arg1` and `arg2` are passed to your script as `ARGS[1]` and `ARGS[2]` of type *String*. Please see [julia docs](https://docs.julialang.org/en/stable/manual/getting-started/) for more options.

(3) To run code in **IJulia notebook**, please see [IJulia](https://github.com/JuliaLang/IJulia.jl).
