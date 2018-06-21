InstallMethod(LiftToSL2Z, [IsMatrix, IsPosInt], function(M, n)
  local F2, S, T, SL2Z, MatS, MatT, MatS_, MatT_, SL2Zn, w;

  if not M in SL(2, Integers mod n) then
    Error("<M> must be a matrix in SL(2,Z/nZ)");
  fi;


  MatS := [[0,-1],[1,0]]; MatT := [[1,1],[0,1]];
  MatS_ := [[ZmodnZObj(0, n), ZmodnZObj(-1,n)], [ZmodnZObj(1,n), ZmodnZObj(0,n)]];
  MatT_ := [[ZmodnZObj(1, n), ZmodnZObj(1,n)], [ZmodnZObj(0,n), ZmodnZObj(1,n)]];

  SL2Z := Group([MatS,MatT]);
  SL2Zn := Group([MatS_, MatT_]);

  w := Factorization(SL2Zn, M);
  F2 := FreeGroup(2);
  w := ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(w));

  return MappedWord(w, [F2.1, F2.2], [MatS, MatT]);
end);

InstallOtherMethod(LiftToSL2Z, [IsMatrixGroup, IsPosInt], function(G, n)
  local gens, gens_Gamma_n;

  if not IsSubgroup(SL(2, Integers mod n), G) then
    Error("<G> must be a subgroup of SL(2,Z/nZ)");
  fi;

  gens := ShallowCopy(GeneratorsOfGroup(G));
  Apply(gens, g -> LiftToSL2Z(g, n));
  gens_Gamma_n := ShallowCopy(GeneratorsOfGroup(PrincipalCongruenceSubgroup(n)));

  return ModularSubgroup(Concatenation(gens, gens_Gamma_n));
end);
