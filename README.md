ModularGroup
====

Algorithms for computing with finite-index subgroups of SL(2,Z)
----


__Usage:__ Either copy this folder into your GAP package directory (usually something like `/opt/gap/gap4r8/pkg`) and load it via
```GAP
LoadPackage("ModularGroup")
```
(this is the recommended way, since it also loads the inline documentation) or use it directly via
```GAP
Read("ModularGroup/lib/ModularSubgroups.gd");
Read("ModularGroup/lib/ModularSubgroups.gi");
```

__Usage example:__ The main object of this package is called `ModularSubgroup`. It represents a finite-index subgroup of _SL(2,Z)_ which is given by the action of the generator matrices _S_ and _T_ on the right cosets. Such an object can be constructed in two different ways:
- either by specifying the action as permutations `s` and `t` and calling the constructor as follows:
```GAP
G := ModularSubgroup(s,t);
```
- or by providing a list of generator matrices:
```GAP
G := ModularSubgroup([ ... ]);
```

The former option is recommended since the latter implicitly computes the coset graph from the generators which might be time-consuming (for more on this, consult the documentation in `doc/manual.pdf`).

Having constructed a modular subgroup, you can apply the various operations this package implements (such as testing if the given group is a congruence subgroup via `IsCongruenceSubgroup(G)`). For more details and a full list of the provided operations, please refer to the documentation in `doc/manual.pdf`. One explicit example is given below:

```GAP
# it is assumed that the package has been loaded as described above
gap> G := ModularSubgroup(
> (1,2,5,3)(4,8,10,9)(6,11,7,12),
> (1,4)(2,6)(3,7)(5,10)(8,12,9,11)
> );
<modular subgroup of index 12>
gap> IsCongruenceSubgroup(G);
true
gap> GeneratorsOfGroup(G);
[ [ [ 1, -2 ], [ 0, 1 ] ], [ [ 1, 0 ], [ 2, 1 ] ] ]
gap> GeneralizedLevel(G);
2
gap> Cusps(G);
[ infinity, 0, 1 ]
```
