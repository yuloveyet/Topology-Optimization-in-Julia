# Topology-Optimization-in-Julia
Julia Codes for Structral Topology Optimization Design


## Codes
`top_oc`

`top_mma`
## Results
![TO design evolution with MMA](./top_mma/res/des_hist.gif)

## Related Julia Pkgs
Start RPL, and
`]` `add ...`

For scientific computing
- LinearAlgebra
- SparseArrays
- Statistics : mean

For images process
- ImageFiltering: imfilter

For modelling and FEM
- Gmsh
- Gridap

For AD
- Zygote

For optimization
- NLopt

## TODO List
- [ ] Sensivity Analasys using Automatic Differentiation
- [ ] Optimization solved with [NLopt](https://github.com/stevengj/nlopt)

## Acknoledgment
TopOpt Group
[top99.m](https://www.topopt.mek.dtu.dk/Apps-and-software/A-99-line-topology-optimization-code-written-in-MATLAB), [top88.m](https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB), [top99neo.m](https://www.topopt.mek.dtu.dk/Apps-and-software/New-99-line-topology-optimization-code-written-in-MATLAB)

Thanks to Prof. [Krister Svanberg](https://people.kth.se/~krille/) with his [Matlab code](http://www.smoptit.se/) for MMA and GCMMA freely available.
