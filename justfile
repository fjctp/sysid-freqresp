set dotenv-load

# Initialize the project for the first time by installing dependencies
init: (repl "-e 'using Pkg; Pkg.instantiate()'")

# Run examples
example-siso1: (repl "example/siso1.jl")
example-siso2: (repl "example/siso2.jl")

# Run a Julia interactive shell
repl command="":
  # Define LD_LIBRARY_PATH due to nix.
  julia -t4 --project=. {{command}}

