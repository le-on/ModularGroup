InstallMethod(ProjectiveModularSubgroup, [IsPerm, IsPerm], function(sp, tp)
  local G, type;

  if not DefinesProjectiveCosetAction(sp, tp) then
    Error("<s> and <t> do not describe the action of the generators S and T on the cosets of a finite-index subgroup of PSL(2,Z)");
  fi;

  type := NewType(FamilyObj(One(SL(2,Integers))),
    IsObject and
    #IsMatrixGroup and
    IsAttributeStoringRep and
    IsComponentObjectRep and
    IsFinitelyGeneratedGroup and
    IsDefaultProjectiveModularSubgroup);

  G := Objectify(type, rec(
    s := sp,
    t := tp
  ));
  return G;
end);

InstallMethod(DefinesProjectiveCosetAction, [IsPerm, IsPerm], function(s, t)
  local index;

  # check relations
  if s^2 <> () or (s*t)^3 <> () then return false; fi;

  # check if the generated subgroup acts transitively
  # this check can be quite costly, so we do it last
  index := Maximum(LargestMovedPoint([s, s^-1, t, t^-1]), 1);
  return IsTransitive(Group(s,t), [1..index]);
end);

InstallMethod(SAction, [IsProjectiveModularSubgroup], function(G)
  return G!.s;
end);

InstallMethod(TAction, [IsProjectiveModularSubgroup], function(G)
  return G!.t;
end);

InstallMethod(Index, "for a projective modular subgroup", [IsProjectiveModularSubgroup], function(G)
  local s, t;
  s := SAction(G);
  t := TAction(G);
  return Maximum(LargestMovedPoint([s, s^-1, t, t^-1]), 1);
end);

InstallMethod(IsCongruenceSubgroup, [IsProjectiveModularSubgroup], function(G)
  local L, R, N, e, m, S, c, d, a, b, l, r, s;
  L := TAction(G);
  R := L*SAction(G)*L;
  N := GeneralizedLevel(G);

  e := 1;
  while RemInt(N, 2) <> 1 do
    e := e * 2;
    N := N / 2;
  od;
  m := N;
  N := e*m;

  if e = 1 then # N is odd
    return (R^2*L^(-(1/2 mod N)))^3 = ();
  elif m = 1 then # N is a power of 2
    S := L^20*R^(1/5 mod N)*L^-4*R^-1;
    return (L*R^-1*L)^-1*S*L*R^-1*L = S^-1 and
           S^-1*R*S = R^25 and
           (S*R^5*L*R^-1*L)^3 = ();
  else
    c := ChineseRem([e, m], [0, 1]) mod N;
    d := ChineseRem([e, m], [1, 0]) mod N;
    a := L^c;
    b := R^c;
    l := L^d;
    r := R^d;
    s := l^20*r^(1/5 mod e)*l^-4*r^-1;
    return a*r*a^-1*r^-1 = () and
           (a*b^-1*a)^4 = () and
           (a*b^-1*a)^2 = (b^-1*a)^3 and
           (a*b^-1*a)^2 = (b^2*a^(-(1/2 mod m)))^3 and
           (l*r^-1*l)^-1*s*l*r^-1*l = s^-1 and
           s^-1*r*s = r^25 and
           (l*r^-1*l)^2 = (s*r^5*l*r^-1*l)^3;
  fi;
end);

InstallMethod(RightCosetRepresentatives, [IsProjectiveModularSubgroup], function(G)
  local s, t, index, H, SC, iso, reps, i, F2, S, T, PSL2Z, epi, epi_psl2z;
  s := SAction(G);
  t := TAction(G);
  index := Index(G);
  H := Group([s,t]);
  SC := StabChain(H);
  reps := [];
  for i in [1..index] do
    reps[i] := InverseRepresentative(SC, i);
  od;

  epi := EpimorphismFromFreeGroup(H);
  Apply(reps, r -> PreImagesRepresentative(epi, r));

  F2 := FreeGroup("S", "T");
  S := F2.1; T := F2.2;
  PSL2Z := F2 / ParseRelators([S,T], "S^2, (S*T)^3");

  F2 := Source(epi);
  epi_psl2z := GroupHomomorphismByImagesNC(F2, PSL2Z, [F2.1, F2.2], [PSL2Z.1, PSL2Z.2]);
  Apply(reps, r -> Image(epi_psl2z, r));

  return reps;
end);

InstallMethod(GeneratorsOfGroup, [IsProjectiveModularSubgroup], function(G)
  local s, t, F2, S, T, PSL2Z, coset_table, H, index;

  s := SAction(G);
  t := TAction(G);
  # Since GAP can only reconstruct generators of a group which is given by a coset graph if the
  # group is a subgroup of a finitely presented group, we need a presentation of PSL(2,Z).
  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  PSL2Z := F2 / ParseRelators([S,T], "S^2, (S*T)^3");

  # In order to compute generators of a subgroup, GAP needs a coset table in terms of generators of the whole
  # group in a specific form. For details, see
  # https://www.gap-system.org/Manuals/doc/ref/chap47.html#X857F239583AFE0B7 and
  # https://www.gap-system.org/Manuals/doc/ref/chap47.html#X7F7F31E47D7F6EF8
  index := Index(G);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(PSL2Z), coset_table);

  return GeneratorsOfGroup(H);
end);

InstallMethod(GeneralizedLevel, [IsProjectiveModularSubgroup], function(G)
  return Order(TAction(G));
end);

InstallMethod(Cusps, [IsProjectiveModularSubgroup], function(G)
  local t, reps, cycles, relevant_reps, c, MatS, MatT, S, T, F2, cusps;

  t := TAction(G);
  reps := RightCosetRepresentatives(G);
  cycles := Orbits(Group(t), [1..Index(G)]);
  relevant_reps := [];
  for c in cycles do
    Add(relevant_reps, reps[c[1]]);
  od;

  MatS := [[0,-1],[1,0]];
  MatT := [[1,1],[0,1]];
  F2 := FreeGroup("S", "T");
  S := F2.1;
  T := F2.2;

  # reverse the words since the action on the cosets is from the right but
  # the action of (P)SL(2,Z) on the extended upper half plane is from the left.
  Apply(relevant_reps, r -> r^-1);
  Apply(relevant_reps, r -> ObjByExtRep(FamilyObj(S), ExtRepOfObj(r)));
  Apply(relevant_reps, r -> MappedWord(r, [S, T], [S^-1, T^-1]));
  
  Apply(relevant_reps, r -> MappedWord(r, [S, T], [MatS, MatT]));

  cusps := List(relevant_reps, r -> MoebiusTransformation(r, infinity));
  return cusps;
end);

InstallMethod(PrintObj, "for projective modular subgroups", [IsProjectiveModularSubgroup], function(G)
  Print("ProjectiveModularSubgroup( ", SAction(G), ", ", TAction(G)," )");
end);

InstallMethod(ViewObj, "for projective modular subgroups", [IsProjectiveModularSubgroup], function(G)
  Print("<projective modular subgroup of index ", Index(G),">");
end);
