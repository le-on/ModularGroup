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

InstallMethod(RightCosetRepresentatives, [IsProjectiveModularSubgroup], function(G)
  local s, t, index, H, SC, iso, reps, i, F2, S, T, PSL2Z, epi, epi_psl2z;
  s := SAction(G);
  t := TAction(G);
  index := Index(G);
  H := Group([s,t]);
  SC := StabChain(H);
  iso := IsomorphismFpGroup(H);
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

InstallMethod(PrintObj, "for projective modular subgroups", [IsProjectiveModularSubgroup], function(G)
  Print("ProjectiveModularSubgroup( ", SAction(G), ", ", TAction(G)," )");
end);

InstallMethod(ViewObj, "for projective modular subgroups", [IsProjectiveModularSubgroup], function(G)
  Print("<projective modular subgroup of index ", Index(G),">");
end);
