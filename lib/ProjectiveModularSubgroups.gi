InstallMethod(ProjectiveModularSubgroup, [IsPerm, IsPerm], function(sp, tp)
  local G, type, tab, index, l1, l2, l3, l4;

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

  index := Maximum(LargestMovedPoint([sp,tp]), 1);
  l1 := ListPerm(sp, index);    l1[1] := l1[1];
  l2 := ListPerm(sp^-1, index); l2[1] := l2[1];
  l3 := ListPerm(tp, index);    l3[1] := l3[1];
  l4 := ListPerm(tp^-1, index); l4[1] := l4[1];
  tab := [l1, l2, l3, l4];
  StandardizeTable(tab);

  G := Objectify(type, rec(
    s := PermList(tab[1]),
    t := PermList(tab[3])
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

InstallMethod(CosetActionOf, [IsMatrix, IsProjectiveModularSubgroup], function(A, G)
  local w, F2;

  if not A in SL(2, Integers) then
    Error("<A> has to be an element of SL(2,Z)");
  fi;

  # By writing A as a word in S and T we can compute its action on the cosets
  # from the actions of S and T.
  w := STDecomposition(A);
  F2 := FreeGroup(2);
  w := ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(w));
  return MappedWord(w, [F2.1, F2.2], [SAction(G), TAction(G)]);
end);

InstallMethod(IsElementOf, [IsMatrix, IsProjectiveModularSubgroup], function(A, G)
  local p;

  if not A in SL(2,Integers) then
    Error("<A> has to be an element of SL(2,Z)");
  fi;

  p := CosetActionOf(A, G);
  return 1^p = 1;
end);
InstallMethod(\in, "for a finite-index subgroup of PSL(2,Z)", [IsMatrix, IsProjectiveModularSubgroup], 10000, function(A, G)
  return IsElementOf(A, G);
end);

InstallMethod(IsSubset, "for two finite-index subgroups of PSL(2,Z)", [IsProjectiveModularSubgroup, IsProjectiveModularSubgroup], function(H, G)
  local gens, result, g, F2, MatS, MatT;
  if Index(H) > Index(G) then return false; fi;
  gens := ShallowCopy(GeneratorsOfGroup(H));
  F2 := FreeGroup("S", "T");
  MatS := [[0,-1],[1,0]];
  MatT := [[1,1], [0,1]];
  Apply(gens, w -> ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(w)));
  Apply(gens, w -> MappedWord(w, [F2.1, F2.2], [MatS, MatT]));
  result := true;
  for g in gens do
    result := result and (g in G);
  od;
  return result;
end);

InstallMethod(\=, "for two finite-index subgroups of PSL(2,Z)", [IsProjectiveModularSubgroup, IsProjectiveModularSubgroup], function(G, H)
  return (Index(G) = Index(H)) and IsSubgroup(G, H);
end);

InstallMethod(Index, "for a projective modular subgroup", [IsProjectiveModularSubgroup], function(G)
  local s, t;
  s := SAction(G);
  t := TAction(G);
  return Maximum(LargestMovedPoint([s, s^-1, t, t^-1]), 1);
end);

InstallMethod(IsCongruence, [IsProjectiveModularSubgroup], function(G)
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
  local s, t, F2, PSL2Z, S, T, H, index, coset_table;
  s := SAction(G);
  t := TAction(G);

  F2 := FreeGroup("S", "T");
  S := F2.1;
  T := F2.2;
  PSL2Z := F2 / ParseRelators([S, T], "S^2, (S*T)^3");

  index := Index(G);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(PSL2Z), coset_table);

  # this is a slicker way to do this without having to explicitly mention a coset
  # table but computing the stabilizer becomes really slow for large (>1000) index
  #hom := GroupHomomorphismByImagesNC(PSL2Z, Group([s,t]), [PSL2Z.1,PSL2Z.2], [s,t]);
  #H := PreImage(hom, Stabilizer(Image(hom), 1));

  return AsList(RightTransversal(PSL2Z, H));
end);

InstallOtherMethod(GeneratorsOfGroup, [IsProjectiveModularSubgroup], function(G)
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

  Apply(relevant_reps, r -> ObjByExtRep(FamilyObj(S), ExtRepOfObj(r)));
  Apply(relevant_reps, r -> MappedWord(r, [S, T], [MatS, MatT]));

  cusps := List(relevant_reps, r -> MoebiusTransformation(r, infinity));
  return cusps;
end);

InstallMethod(CuspWidth, [IsRat, IsProjectiveModularSubgroup], function(c, G)
  local p, q, gcd, g, w, F2, S, T, PSL2Z, reps, r, s, t, k;

  p := NumeratorRat(c);
  q := DenominatorRat(c);
  gcd := Gcdex(q, -p);
  g := [[gcd.coeff1, p], [gcd.coeff2, q]] * [[0,-1], [1,0]];
  w := STDecomposition(g);

  F2 := FreeGroup("S", "T");
  S := F2.1; T := F2.2;
  PSL2Z := F2 / ParseRelators([S,T], "S^2, (S*T)^3");
  w := ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(w)));

  reps := ShallowCopy(RightCosetRepresentatives(G));
  Apply(reps, r -> ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(r))));
  s := SAction(G);
  t := TAction(G);
  k := 1;
  for r in reps do
    if 1^MappedWord(w*r^-1, [PSL2Z.1,PSL2Z.2], [s,t]) = 1 then
      break;
    fi;
    k := k+1;
  od;

  return CycleLength(t, [1..Index(G)], k);
end);
InstallOtherMethod(CuspWidth, [IsInfinity, IsProjectiveModularSubgroup], function(c, G)
  return CycleLength(TAction(G), [1..Index(G)], 1);
end);

InstallMethod(CuspsEquivalent, [IsRat, IsRat, IsProjectiveModularSubgroup], function(c1, c2, G)
  local p1, p2, q1, q2, gcd1, gcd2, g1, g2, w1, w2, F2, S, T, PSL2Z, reps, r, s, t, k1, k2;

  # find a matrix g1 that maps infinity to c1 and express it as a word w1 in S and T
  p1 := NumeratorRat(c1);
  q1 := DenominatorRat(c1);
  gcd1 := Gcdex(q1, -p1);
  g1 := [[gcd1.coeff1, p1], [gcd1.coeff2, q1]] * [[0,-1], [1,0]];
  w1 := STDecomposition(g1);

  # find a matrix g2 that maps infinity to c2 and express it as a word w2 in S and T
  p2 := NumeratorRat(c2);
  q2 := DenominatorRat(c2);
  gcd2 := Gcdex(q2, -p2);
  g2 := [[gcd2.coeff1, p2], [gcd2.coeff2, q2]] * [[0,-1], [1,0]];
  w2 := STDecomposition(g2);

  F2 := FreeGroup("S", "T");
  S := F2.1; T := F2.2;
  PSL2Z := F2 / ParseRelators([S,T], "S^2, (S*T)^3");
  w1 := ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(w1)));
  w2 := ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(w2)));

  reps := ShallowCopy(RightCosetRepresentatives(G));
  Apply(reps, r -> ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(r))));
  s := SAction(G);
  t := TAction(G);

  # find the coset in which w1 lies
  k1 := 1;
  for r in reps do
    if 1^MappedWord(w1*r^-1, [PSL2Z.1,PSL2Z.2], [s,t]) = 1 then
      break;
    fi;
    k1 := k1+1;
  od;

  # find the coset in which w2 lies
  k2 := 1;
  for r in reps do
    if 1^MappedWord(w2*r^-1, [PSL2Z.1,PSL2Z.2], [s,t]) = 1 then
      break;
    fi;
    k2 := k2+1;
  od;

  # c1 and c2 are equivalent modulo G if and only if the cosets in which w1 and w2 lie
  # are in the same cycle of the T-action
  return k2 in Orbit(Group(t), k1);
end);
InstallOtherMethod(CuspsEquivalent, [IsRat, IsInfinity, IsProjectiveModularSubgroup], function(c1, c2, G)
  local p, q, gcd, g, w, F2, S, T, PSL2Z, reps, r, s, t, k;

  p := NumeratorRat(c1);
  q := DenominatorRat(c1);
  gcd := Gcdex(q, -p);
  g := [[gcd.coeff1, p], [gcd.coeff2, q]] * [[0,-1], [1,0]];
  w := STDecomposition(g);

  F2 := FreeGroup("S", "T");
  S := F2.1; T := F2.2;
  PSL2Z := F2 / ParseRelators([S,T], "S^2, (S*T)^3");
  w := ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(w)));

  reps := ShallowCopy(RightCosetRepresentatives(G));
  Apply(reps, r -> ElementOfFpGroup(FamilyObj(One(PSL2Z)), ObjByExtRep(FamilyObj(S), ExtRepOfObj(r))));
  s := SAction(G);
  t := TAction(G);
  k := 1;
  for r in reps do
    if 1^MappedWord(w*r^-1, [PSL2Z.1,PSL2Z.2], [s,t]) = 1 then
      break;
    fi;
    k := k+1;
  od;

  return k in Orbit(Group(t), 1);
end);
InstallOtherMethod(CuspsEquivalent, [IsInfinity, IsRat, IsProjectiveModularSubgroup], function(c1, c2, G)
  return CuspsEquivalent(c2, c1, G);
end);
InstallOtherMethod(CuspsEquivalent, [IsInfinity, IsInfinity, IsProjectiveModularSubgroup], function(c1, c2, G)
  return true;
end);

InstallMethod(CosetRepresentativeOfCusp, [IsRat, IsProjectiveModularSubgroup], function(c, G)
  local t, reps, cycles, relevant_reps, d, MatS, MatT, S, T, F2, r;

  t := TAction(G);
  reps := RightCosetRepresentatives(G);
  cycles := Orbits(Group(t), [1..Index(G)]);
  relevant_reps := [];
  for d in cycles do
    Add(relevant_reps, reps[d[1]]);
  od;

  MatS := [[0,-1],[1,0]];
  MatT := [[1,1],[0,1]];
  F2 := FreeGroup("S", "T");
  S := F2.1;
  T := F2.2;

  Apply(relevant_reps, r -> ObjByExtRep(FamilyObj(S), ExtRepOfObj(r)));

  for r in relevant_reps do
    if CuspsEquivalent(c, MoebiusTransformation(MappedWord(r, [S, T], [MatS, MatT]), infinity), G) then
      return r;
    fi;
  od;
end);
InstallOtherMethod(CosetRepresentativeOfCusp, [IsInfinity, IsModularSubgroup], function(c, G)
  return One(FreeGroup(["S", "T"]));
end);

InstallMethod(LiftToSL2ZEven, [IsProjectiveModularSubgroup], function(G)
  return ModularSubgroup(SAction(G), TAction(G));
end);

InstallMethod(LiftToSL2ZOdd, [IsProjectiveModularSubgroup], function(G)
  local gens, F2, S, T;
  gens := ShallowCopy(GeneratorsOfGroup(G));
  F2 := FreeGroup("S", "T");
  S := [[0,-1],[1,0]];
  T := [[1,1],[0,1]];
  Apply(gens, g -> MappedWord(ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(g)), [F2.1, F2.2], [S, T]));
  return ModularSubgroup(gens);
end);

InstallMethod(IndexModN, [IsProjectiveModularSubgroup, IsPosInt], function(G, N)
  local lift;
  lift := LiftToSL2ZEven(G);
  return IndexModN(lift, N);
end);

InstallMethod(Deficiency, [IsProjectiveModularSubgroup], function(G)
  if IsCongruence(G) then return 1; fi;
  return Index(G) / IndexModN(G, GeneralizedLevel(G));
end);

InstallMethod(Deficiency, [IsProjectiveModularSubgroup, IsPosInt], function(G, N)
  return Index(G) / IndexModN(G, N);
end);

InstallMethod(Conjugate, [IsProjectiveModularSubgroup, IsMatrix], function(G, A)
  local a, s, t;
  a := CosetActionOf(A, G);
  s := SAction(G);
  t := TAction(G);
  return ProjectiveModularSubgroup(a^-1*s*a, a^-1*t*a); # this is A^-1*G*A
  # note that if G = (\sigma_S, \sigma_T) then A^-1*G*A = (\sigma_A*\sigma_S*\sigma_A^-1, \sigma_A*\sigma_T*\sigma_A^-1)
  # but GAP multiplies permutations the other way around
end);

InstallMethod(NormalCore, [IsProjectiveModularSubgroup], function(G)
  local s, t, F2, S, T, PSL2Z, index, m, Sd, core, coset_table;

  s := SAction(G);
  t := TAction(G);

  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  PSL2Z := F2 / ParseRelators([S,T], "S^2, (S*T)^3");
  S := PSL2Z.1;
  T := PSL2Z.2;

  index := Index(G);
  Sd := SymmetricGroup(index);

  m := GroupHomomorphismByImagesNC(PSL2Z, Sd, [S,T], [s,t]);
  core := Kernel(m);
  coset_table := CosetTableBySubgroup(PSL2Z, core);
  return ProjectiveModularSubgroup(PermList(coset_table[1]), PermList(coset_table[3]));
end);

InstallMethod(QuotientByNormalCore, [IsProjectiveModularSubgroup], function(G)
  local NC, s, t, F2, S, T, PSL2Z, index, coset_table, H;

  NC := NormalCore(G);
  s := SAction(NC);
  t := TAction(NC);

  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  PSL2Z := F2 / ParseRelators([S,T], "S^4, (S*T)^3");

  index := Index(NC);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(PSL2Z), coset_table);

  return FactorGroupNC(PSL2Z, H);
end);

InstallMethod(AssociatedCharacterTable, [IsProjectiveModularSubgroup], function(G)
  local F;
  F := QuotientByNormalCore(G);
  return CharacterTable(F);
end);

InstallMethod(Genus, [IsProjectiveModularSubgroup], function(G)
  local reps, F2, S, T, t, e, v, glued, i, j, r1, r2, c;
  reps := ShallowCopy(RightCosetRepresentatives(G));
  F2 := FreeGroup(["S", "T"]);
  S := [[0,-1], [1,0]]; T := [[1,1], [0,1]];
  Apply(reps, w -> ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(w)));
  Apply(reps, w -> MappedWord(w, [F2.1, F2.2], [S, T]));

  t := Length(reps);
  glued := ListWithIdenticalEntries(t, false);
  c := 0;
  for i in [1..t] do
    if glued[i] then continue; fi;
    r1 := reps[i];
    for j in [1..t] do
      if glued[j] then continue; fi;
      r2 := reps[j];
      if IsElementOf(r1*S*r2^-1, G) then
        glued[j] := true;
        if i = j then c := c+1; fi;
        break;
      fi;
    od;
    glued[i] := true;
  od;
  t := t+c;
  e := 3*t/2;
  v := Length(Cycles(TAction(G), [1..Index(G)])) # number of cusps
     + Length(Cycles(TAction(G)^-1*SAction(G), [1..Index(G)])) # number of cycles in T^-1*S
     + Length(Filtered(Cycles(SAction(G), [1..Index(G)]), l -> Length(l) = 1)); # number of singleton cycles of S

  return (2-(v-e+t))/2;
end);

InstallMethod(IO_Pickle, "for a projective modular subgroup", [IsFile, IsProjectiveModularSubgroup], function(f, G)
  local nr;
  nr := IO_AddToPickled(G); # this should always be false (?)
  if IO_Write(f, "PMOD") = fail then IO_FinalizePickled(); return IO_Error; fi;
  if IO_Pickle(f, SAction(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  if IO_Pickle(f, TAction(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;

  if HasIsCongruence(G) then
    if IO_Pickle(f, IsCongruence(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  if HasCusps(G) then
    if IO_Pickle(f, Cusps(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  ## There is no implementation of IO_Pickle for words in a finitely presented group
  #if HasRightCosetRepresentatives(G) then
  #  if IO_Pickle(f, RightCosetRepresentatives(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  #else
  #  if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  #fi;

  if HasGeneralizedLevel(G) then
    if IO_Pickle(f, GeneralizedLevel(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  ## There is no implementation of IO_Pickle for words in a finitely presented group
  #if HasGeneratorsOfGroup(G) then
  #  if IO_Pickle(f, GeneratorsOfGroup(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  #else
  #  if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  #fi;

  if HasNormalCore(G) then
    if IO_Pickle(f, NormalCore(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  if HasQuotientByNormalCore(G) then
    if IO_Pickle(f, QuotientByNormalCore(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  if HasAssociatedCharacterTable(G) then
    if IO_Pickle(f, AssociatedCharacterTable(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  if HasGenus(G) then
    if IO_Pickle(f, Genus(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  IO_FinalizePickled();
  return IO_OK;
end);

IO_Unpicklers.PMOD := function(f)
  local H, s, t, cong, cusps, level, normal_core, quotient_nc, ctbl, genus;
  s := IO_Unpickle(f);
  if s = IO_Error then return IO_Error; fi;
  t := IO_Unpickle(f);
  if t = IO_Error then return IO_Error; fi;
  H := ProjectiveModularSubgroup(s, t);
  IO_AddToUnpickled(H);

  cong := IO_Unpickle(f);
  if cong = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  cusps := IO_Unpickle(f);
  if cusps = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  level := IO_Unpickle(f);
  if level = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  normal_core := IO_Unpickle(f);
  if normal_core = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  quotient_nc := IO_Unpickle(f);
  if quotient_nc = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  ctbl := IO_Unpickle(f);
  if ctbl = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  genus := IO_Unpickle(f);
  if genus = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;

  if cong <> fail then SetIsCongruence(H, cong); fi;
  if cusps <> fail then SetCusps(H, cusps); fi;
  if level <> fail then SetGeneralizedLevel(H, level); fi;
  if normal_core <> fail then SetNormalCore(H, normal_core); fi;
  if quotient_nc <> fail then SetQuotientByNormalCore(H, quotient_nc); fi;
  if ctbl <> fail then SetAssociatedCharacterTable(H, ctbl); fi;
  if genus <> fail then SetGenus(H, genus); fi;

  IO_FinalizeUnpickled();
  return H;
end;

InstallMethod(PrintObj, "for projective modular subgroups", [IsProjectiveModularSubgroup], function(G)
  Print("ProjectiveModularSubgroup( ", SAction(G), ", ", TAction(G)," )");
end);

InstallMethod(ViewObj, "for projective modular subgroups", [IsProjectiveModularSubgroup], function(G)
  Print("<projective modular subgroup of index ", Index(G),">");
end);
