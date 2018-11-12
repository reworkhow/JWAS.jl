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

To use the latest/beta features under development, run `Pkg.add(PackageSpec(name="JWAS", rev="master"))` to get the
newest unofficial JWAS. Run `Pkg.free("JWAS")` to go back to the official one.

#### Jupyter Notebook

If you prefer “reproducible research”, an interactive Jupyter Notebook interface is available
for Julia (and therefore JWAS). The Jupyter Notebook is an open-source web application for creating
and sharing documents that contain live code, equations, visualizations and explanatory text.
To install IJulia for Jupyter Notebook, please go to [IJulia](https://github.com/JuliaLang/IJulia.jl).

#### Docker


!!! note "Jupyter Notebooks with JWAS via Docker"

    Docker provides a straightforward way to install Jupyter Notebooks with JWAS.


- Install Docker from [here](https://docs.docker.com/install/) for your platform.

- From a terminal (on Mac or Linux), run the command:

```bash
docker run -it --rm -p 8888:8888 qtlrocks/jwas-docker
```

This will start a Jupyter Notebook server listening for HTTP connections on port 8888 with a randomly generated authentication
token. Examples for JWAS can be accessed from the notebook: `notebooks/0_index.ipynb`.

The directories and files created within the Docker container will be lost when the container is stopped. To save your work
on the host machine, a directory on the host machine can be mounted as a folder in the container with the `-v` option. After `cd` into your
working directory on your local machine or a server, run the command

```bash
docker run -it --rm -p 8888:8888 -v `pwd`:/home/jovyan/work qtlrocks/jwas-docker
```

This command creates a Docker container with the folder `/home/jovyan/work` with the contents of `pwd` of the host machine. Files and
directories that are in the folder `pwd` will not be lost when the container is stopped.  

After running this command, it is expected to prompt something like

```
[I 10:41:54.774 NotebookApp] Writing notebook server cookie secret to /home/ubuntu/.local/share/jupyter/runtime/notebook_cookie_secret
[I 10:41:54.920 NotebookApp] Serving notebooks from local directory: /home/ubuntu
[I 10:41:54.920 NotebookApp] 0 active kernels
[I 10:41:54.920 NotebookApp] The Jupyter Notebook is running at:
[I 10:41:54.920 NotebookApp] http://0.0.0.0:8888/?token=75ad671f75b4c47be70591f46bec604997d8a9bd9dd51f0d
[I 10:41:54.920 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 10:41:54.921 NotebookApp]

    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://0.0.0.0:8888/?token=75ad671f75b4c47be70591f46bec604997d8a9bd9dd51f0d
```

Then, open the url in an internet browser (IE, Firefox, Chrome, Safari, etc) if JWAS-docker is launched on your local machine.

If you prefer running scripts using linux commands in Bash instead of Jupyter Notebook, please run the command

```bash
docker run -it --rm -v `pwd`:/home/jovyan/work qtlrocks/jwas-docker bash
```


#### Standalone application

!!! note "standalone application (no installation required)"

    A fully self-contained application for JWAS (no installation required) will come out next year.


## Access documentation

!!! warning

    Please load the JWAS package at first.

To show the basic information (README file) of JWAS in REPL or IJulia notebook using `?JWAS`
and press enter.

For help on a specific function, type ? followed by its name, e.g. `?runMCMC` and press enter
in REPL or IJulia notebook.

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

(3) To run code in **Jupyter Notebook**, please see [IJulia](https://github.com/JuliaLang/IJulia.jl).

(4) To run code in **Jupyter Notebook via Docker**, please see [Docker](@ref).
