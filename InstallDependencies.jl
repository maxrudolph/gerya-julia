# Install dependencies
# Installs the packages listed in Project.toml into this repository's environment.
# Run from anywhere:  julia InstallDependencies.jl
# Then run scripts with the same environment, e.g.  julia --project PlumeModel_lower_eta.jl ...
# (Jupyter kernels started with --project=@. in this directory use it automatically.)
#
# Pardiso (MKL) is not included: a few notebooks (M04-YinYang, M06-Plume, Chapter13-I2ELVIS) use it,
# but MKL is not available on every platform (e.g. Apple silicon). To use it where MKL is available:
#     julia --project -e 'using Pkg; Pkg.add("Pardiso")'
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
