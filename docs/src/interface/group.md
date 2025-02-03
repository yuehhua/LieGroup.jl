# An interface for Lie groups

```@docs
LieGroup
AbstractLieAlgebraTangentVector
AbstractLieGroupPoint
```

## Functions on Lie groups

```@autodocs
Modules = [LieGroups]
Pages = ["src/interface.jl"]
Order = [:function]
```

## Internal functions and macros

```@docs
LieGroups.CommonUnitarySubAlgebras
LieGroups.@default_lie_group_fallbacks
```

## Literature

```@bibliography
Pages = ["group.md"]
Canonical=false
```