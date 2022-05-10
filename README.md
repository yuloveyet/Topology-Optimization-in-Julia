# Topology-Optimization-in-Julia

Julia Codes for Structural Topology Optimization Design

## Codes

A personal [Julia](https://epubs.siam.org/doi/10.1137/141000671) code is given mainly based on a compact and efficient Matlab implementation of compliance topology optimization (TO) for 2D continua[^1], which is a v3.0 version of the celebrated $\text{top99}$ Matlab code developed by Sigmund[^2] and $\text{top88}$ by its heir[^3].

Assemble just one half of the sysmetric stiffness matrix, thus substantial speedups are acheived.

`top_oc` and `top_mma` contain corresponding files related to the TO with OC and MMA algorithms, respectively.

Running codes could be tried out as:

```julia
include("./top99neo_mma.jl")
setup = SetUp() # problem setup
mat = Mat() # material property
disfeature = DiscretizationFeature(setup, mat) # model discretization
load = LoadsSupportsBCs(setup, disfeature) # boudary conditions
ini = Initialization(setup, disfeature, mat) # initial conditions
filter = Filter(setup) # filtering
xPhys, opt_hist, vf_hist, anim = Optimization(setup, mat, load, filter, ini, disfeature) # optimization process
gif(anim, "./res/des_hist.gif", fps=20) # design result visulization
```

⚠️ Disclaimer: The authors reserves all rights but do not guaranty that the code is free from errors. Furthermore, we shall not be liable in any event.

## Results

A benchmark MBB example is presented. TO design result evolution is shown as below.

<!-- ![TO design evolution with MMA](./top_mma/res/des_hist.gif) -->

<div align=center>
<img src=./top_mma/res/des_hist.gif width="500">
👍 💯
</div>

## Required Julia Packages

Run Julia REPL, enter `]` to bring up Julia's [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/),
and add the listed packages:

> julia> ]
> (@v1.7) pkg> add #pkg_name#

- Scientific computing
  - LinearAlgebra
  - SparseArrays
  - Statistics : mean
- Image process
  - ImageFiltering: imfilter
- Modelling
  - [Gmsh](https://onlinelibrary.wiley.com/doi/10.1002/nme.2579)
- FEM
  - [Gridap](https://joss.theoj.org/papers/10.21105/joss.02520)
- AD
  - [ForwardDiff](https://arxiv.org/abs/1607.07892)
  - [Zygote](https://arxiv.org/abs/1810.07951)
- Optimization
  - [NLopt](https://github.com/stevengj/nlopt)
- Visualization
  - Plots

## TODO List

- [x] [`top99neo.jl`](./top_oc/top99neo.jl)$\text{top99neo.m}$ rewritten in Julia
- [x] [`MMA.jl`](./top_mma/MMA.jl)MMA algorithm ($\text{mmasub.m + subsolve.m}$) rewritten in Julia
- [x] [`top99neo_mma`](./top_mma/top99neo_mma.jl)$\text{top99neo +}$ MMA in Julia
- [ ] `top99neo_AD`Sensitivity Analysis using Automatic Differentiation
- [ ] `top99neo_NLopt`Optimization solved with NLopt

## Acknowledgements

TopOpt Group 🇩🇰
Matlab codes for topology optimization

- v1.0 [top99.m](https://www.topopt.mek.dtu.dk/Apps-and-software/A-99-line-topology-optimization-code-written-in-MATLAB)
    Educatianal TO enlightenment for every beginners
- v2.0 [top88.m](https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB)
    Loop vectorization and memory preallocation
- v3.0 [top99neo.m](https://www.topopt.mek.dtu.dk/Apps-and-software/New-99-line-topology-optimization-code-written-in-MATLAB)
   Half matrix assembly operation, filter implementation and volume-preserving density projection

Prof. [Krister Svanberg](https://people.kth.se/~krille/) 🇸🇪
freely available [Matlab code](http://www.smoptit.se/) for CCSA/MMA[^4] and GCMMA

## Author ©️

Yu Li 📧 liyu_npu@outlook.com

🇨🇳 🎓 1️⃣9️⃣9️⃣0️⃣  ♑

Hobby 🎧 🃏 🎮 🏀 🏊 🏃 🚴‍♂️

Food 🍦 🦞 🍣 🌽 🍌

**References**
[^1]:Sigmund, O. (2001). A 99 line topology optimization code written in Matlab. Structural and multidisciplinary optimization, 21(2), 120-127.
[^2]: Ferrari, F., & Sigmund, O. (2020). A new generation 99 line Matlab code for compliance topology optimization and its extension to 3D. Structural and Multidisciplinary Optimization, 62(4), 2211-2228.
[^3]:Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B. S., & Sigmund, O. (2011). Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization, 43(1), 1-16.
[^4]: Svanberg, K. (2002). A class of globally convergent optimization methods based on conservative convex separable approximations. SIAM journal on optimization, 12(2), 555-573.
