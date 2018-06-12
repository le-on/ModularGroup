gap> G := ProjectiveModularSubgroup((1,2)(3,4)(5,6)(7,8)(9,10), (1,4)(2,5,9,10,8)(3,7,6));
<projective modular subgroup of index 10>
gap> IsCongruenceSubgroup(G);
false
gap> Cusps(G);
[ infinity, 0, 1 ]
gap> CuspWidth(infinity, G);
2
gap> CuspWidth(1, G);
3
gap> CuspWidth(0, G);
5
gap> GeneralizedLevel(G);
30
gap> Deficiency(G, 30);
10
gap> RightCosetRepresentatives(G);
[ <identity ...>, S, T, S*T, S*T^-1, T*S, S*T*S, S*T^2, S*T^-1*S, S*T^-2 ]
