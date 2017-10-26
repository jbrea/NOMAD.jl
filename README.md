# The NOMAD module for Julia
This module provides a [Julia-language](http://julialang.org/) interface to the 
free/open-source [NOMAD](https://www.gerad.ca/nomad/Project/Home.html) library 
for nonlinear optimization of black-box functions. Currently only a limited set
of features of [NOMAD](https://www.gerad.ca/nomad/Project/Home.html) is
available through this interface.

## Installation

Currently the installation is only tested on Linux. Please adapt
[deps/build.jl](deps/build.jl) for other platforms.

Within Julia, you can install the NOMAD.jl package with the package manager: 
`Pkg.clone("https://github.com/jbrea/NOMAD.jl")`

## Basic Usage
```julia
using NOMAD
function rosenbrock(x)
    res = 0.
    for i in 1:length(x) - 1
		res += 100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2
	end
	res
end
ev = Evaluator(rosenbrock, [], 5)
opt = Opt(ev, UPPER_BOUND = 10*ones(5), LOWER_BOUND = -10*ones(5), 
		  X0 = zeros(5), DISPLAY_DEGREE = 0)
optimize(opt)
```

```julia
using NOMAD
ev = Evaluator((x, y) -> - x^2 + (y - 1)^2, [(x, y) -> x - y], 2, vectorin = false)
opt = Opt(ev, UPPER_BOUND = [5, 6], LOWER_BOUND = [-2, -6], X0 = [1, 2])
optimize(opt)
```

```julia
using NOMAD
ev = Evaluator(x -> x[5], [(:PB, x -> (norm(x-1)^2 - 25)), (:EB, x -> 25 - norm(x + 1)^2
)], 5)
opt = Opt(ev, UPPER_BOUND = [5, 6, 7, Inf64, Inf64], LOWER_BOUND = -6*ones(5), X0 = zero
s(5), MAX_BB_EVAL = 100, SEED = -1)
optimize(opt)
```

## Using JuMP
NOMAD implements the [MathProgBase
interface](http://mathprogbasejl.readthedocs.org/en/latest/nlp.html) for
nonlinear optimization, which means that it can be used interchangeably with
other optimization packages from modeling packages like
[JuMP](https://github.com/JuliaOpt/JuMP.jl).

```julia
using JuMP, NOMAD
m = Model(solver = NOMAD.NOMADSolver())
a1 = 2
b1 = 0
a2 = -1
b2 = 1

@variable(m, x1)
@variable(m, x2 >= 0)

@NLobjective(m, Min, sqrt(x2))
@NLconstraint(m, x2 >= (a1*x1+b1)^3)
@NLconstraint(m, x2 >= (a2*x1+b2)^3)

setvalue(x1, 1.234)
setvalue(x2, 5.678)

status = solve(m)
println("got ", getobjectivevalue(m), " at ", [getvalue(x1),getvalue(x2)])
```

## Help on NOMAD Parameters		

```julia
NOMAD.help()
NOMAD.help(:UPPER_BOUND)
```
