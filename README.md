
# PSA.jl

## About

PSA.jl is a partial implementation of [Python for Power System
Analysis (PyPSA)](https://github.com/PyPSA/PyPSA) in
the programming language [Julia](https://julialang.org/).

PSA.jl has been created primarily to take advantage of the speed,
readability and features of the optimisation framework
[JuMP](https://github.com/JuliaOpt/JuMP.jl).

PSA.jl does not yet exist independently of PyPSA, in that you have to
build your network first in PyPSA and export it in CSV format before
PSA.jl can import the network and work with it.

So far the following functionality is implemented in PSA.jl:

* Import from CSV folder
* Export to CSV folder
* Most features of Linear Optimal Power Flow (LOPF), i.e. total electricity/energy system least-cost investment optimisation

What is not yet implemented, but coming soon:

* Normal power flow
* Unit commitment
* Stores
* LOPF with non-unity snapshot weightings
* Links with more than one output
* Standing losses in storage
* Some time-dependent efficiencies and max/min limits
* Import/Export from/to NetCDF


## Required packages

* JuMP.jl
* CSV.jl
* LightGraphs.jl
* DataFrames.jl
* PyCall.jl
* Gurobi.jl

## Basic usage

To install PSA.jl, execute `git clone git@github.com:PyPSA/PSA.jl.git`
in the folder `/path/to/lib` of your choice.

Then in your Julia script:

```
push!(LOAD_PATH, "/path/to/lib")

import PSA

network = PSA.import_network("/your/exported/PyPSA/network/")

using Clp

solver = ClpSolver()

m = PSA.lopf(network, solver)

print(m.objVal)
```


## Licence

Copyright 2017-2018 Fabian Hofmann (FIAS), Tom Brown (KIT IAI)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either [version 3 of the
License](LICENSE.txt), or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
[GNU General Public License](LICENSE.txt) for more details.
