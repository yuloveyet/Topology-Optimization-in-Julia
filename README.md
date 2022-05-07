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

For computing
- LinearAlgebra
- SparseArrays
- Statistics : mean

For images process
- ImageFiltering: imfilter

## TODO List
- [ ] Sensivity Analasys using Automatic Differentiation
- [ ] Optimization solved with [NLopt](https://github.com/stevengj/nlopt)

## Acknoledgment
TopOpt Group
top99, top88, [top99neo](https://www.topopt.mek.dtu.dk/Apps-and-software/New-99-line-topology-optimization-code-written-in-MATLAB)

Thanks to Prof. [Krister Svanberg](https://people.kth.se/~krille/) with his [Matlab code](http://www.smoptit.se/) for MMA and GCMMA freely available.
