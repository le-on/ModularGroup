gap> G := ModularSubgroup([
> [[1,2],[0,1]],
> [[1,0],[2,1]]
> ]);
<modular subgroup of index 12>
gap> IsCongruenceSubgroup(G);
true
gap> Cusps(G);
[ infinity, 0, 1 ]
gap> RightCosetRepresentatives(G);
[ <identity ...>, S, S^-1, T, S^2, S*T, S^-1*T, T*S, T*S^-1, S^2*T, S*T*S,
  S*T*S^-1 ]
gap> GeneralizedLevel(G);
2
gap> GeneratorsOfGroup(G);
[ T^-2, S*T^-2*S^-1 ]
gap> MatrixGeneratorsOfGroup(G);
[ [ [ 1, 2 ], [ 0, 1 ] ], [ [ 1, 0 ], [ 2, 1 ] ] ]
gap> SAction(G);
(1,2,5,3)(4,8,10,9)(6,11,7,12)
gap> TAction(G);
(1,4)(2,6)(3,7)(5,10)(8,12,9,11)
gap> CuspWidth(infinity, G);
2
gap> CuspWidth(1, G);
2
gap> CuspWidth(0, G);
2
gap> IndexModN(G, 4);
12
gap> IndexModN(G, 2);
6
gap> Deficiency(G, 2);
2
gap> Deficiency(G, 4);
1
gap> Projection(G);
<projective modular subgroup of index 6>
