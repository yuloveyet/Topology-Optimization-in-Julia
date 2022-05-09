# Topology-Optimization-in-Julia

[Julia](https://epubs.siam.org/doi/10.1137/141000671) Codes for Structural Topology Optimization Design

## Codes

`top_oc`

`top_mma`

## Results

![TO design evolution with MMA](./top_mma/res/des_hist.gif)

👍 💯

## Required Julia Packages

Run Julia REPL, enter `]` to bring up Julia's [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/),
and add the listed packages:

```julia
julia> ]
(@v1.7) pkg> add ...
```

For scientific computing

- LinearAlgebra
- SparseArrays
- Statistics : mean

For images process

- ImageFiltering: imfilter

For modelling and FEM

- Gmsh
- [Gridap](https://joss.theoj.org/papers/10.21105/joss.02520)

For AD

- [ForwardDiff](https://arxiv.org/abs/1607.07892)
- [Zygote](https://arxiv.org/abs/1810.07951)

For optimization

- NLopt

For visualization

- Plots

## TODO List

- [x] top99neo.m rewritten in Julia
- [x] MMA rewritten in Julia
- [x] top_mma in Julia
- [ ] Sensitivity Analysis using Automatic Differentiation
- [ ] Optimization solved with [NLopt](https://github.com/stevengj/nlopt)

## Acknowledgements

TopOpt Group 🇩🇰
[top99.m](https://www.topopt.mek.dtu.dk/Apps-and-software/A-99-line-topology-optimization-code-written-in-MATLAB), [top88.m](https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB), [top99neo.m](https://www.topopt.mek.dtu.dk/Apps-and-software/New-99-line-topology-optimization-code-written-in-MATLAB)

Thanks to Prof. [Krister Svanberg](https://people.kth.se/~krille/) 🇸🇪 with his [Matlab code](http://www.smoptit.se/) for [CCSA/MMA][1] and GCMMA freely available.

## References

[1]: Svanberg, K. (2002). A class of globally convergent optimization methods based on conservative convex separable approximations. SIAM journal on optimization, 12(2), 555-573.

## Author ©️

🎓 Yu Li 👓
🇨🇳 

♑

1️⃣9️⃣9️⃣0️⃣ 🐴

Hobby 🎧 🃏 🎮 🏀 🏊 🏃 🚴‍♂️

Food 🍦 🦞 🍣 🌽 🍌



📧 liyu_npu@outlook.com


