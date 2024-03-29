InstallMethod(ModularSubgroup, [IsPerm, IsPerm], function(sp, tp)
  local G, type, tab, index, l1, l2, l3, l4;

  if not DefinesCosetActionST(sp, tp) then
    Error("<s> and <t> do not describe the action of the generators S and T on the cosets of a finite-index subgroup of SL(2,Z)");
  fi;

  type := NewType(FamilyObj(One(SL(2,Integers))),
    IsObject and
    IsMatrixGroup and
    IsAttributeStoringRep and
    IsComponentObjectRep and
    IsFinitelyGeneratedGroup and
    IsDefaultModularSubgroup);

  index := Maximum(LargestMovedPoint([sp,tp]), 1);
  l1 := ListPerm(sp, index);    l1[1] := l1[1];
  l2 := ListPerm(sp^-1, index); l2[1] := l2[1];
  l3 := ListPerm(tp, index);    l3[1] := l3[1];
  l4 := ListPerm(tp^-1, index); l4[1] := l4[1];
  tab := [l1, l2, l3, l4];
  StandardizeTable(tab);

  sp := PermList(tab[1]);
  tp := PermList(tab[3]);

  G := Objectify(type, rec(
    s := sp,
    t := tp,
    r := sp^-1*tp^-1*sp,
    j := sp^-1*tp^-1
  ));
  return G;
end);

InstallMethod(ModularSubgroupST, [IsPerm, IsPerm], function(sp, tp)
  return ModularSubgroup(sp, tp);
end);

InstallMethod(ModularSubgroupRT, [IsPerm, IsPerm], function(rp, tp)
  local sp;
  sp := rp*tp^-1*rp;
  return ModularSubgroup(sp, tp);
end);

InstallMethod(ModularSubgroupSJ, [IsPerm, IsPerm], function(sp, jp)
  local tp;
  tp := jp^-1*sp^-1;
  return ModularSubgroup(sp, tp);
end);

InstallOtherMethod(ModularSubgroup, [IsList], function(gens)
  local G, a;
  if Length(gens) = 0 then
    Error("<gens> must not be empty");
  fi;
  if not (IsMatrix(gens[1]) and gens[1] in SL(2,Integers)) then
    Error("<gens> needs to be a non-empty list of matrices in SL(2,Z)");
  fi;

  a := CosetActionFromGenerators(gens);
  G := ModularSubgroup(a[1], a[2]);
  SetMatrixGeneratorsOfGroup(G, gens);
  return G;
end);

InstallMethod(DefinesCosetActionST, [IsPerm, IsPerm], function(s, t)
  local index;

  # check relations
  if s^4 <> () or (s^3*t)^3 <> () or s^2*t*s^-2*t^-1 <> () then return false; fi;

  # check if the generated subgroup acts transitively
  index := Maximum(LargestMovedPoint([s, t]), 1);
  return IsTransitive(Group(s,t), [1..index]);
end);

InstallMethod(DefinesCosetActionRT, [IsPerm, IsPerm], function(r, t)
  local index;

  # check relations
  if (r*t^-1*r)^4 <> () or ((r*t^-1*r)^3*t)^3 <> () or (r*t^-1*r)^2*t*(r*t^-1*r)^-2*t^-1 <> () then return false; fi;

  # check if the generated subgroup acts transitively
  index := Maximum(LargestMovedPoint([r, t]), 1);
  return IsTransitive(Group(r,t), [1..index]);
end);

InstallMethod(DefinesCosetActionSJ, [IsPerm, IsPerm], function(s, j)
  local index;

  # check relations
  if s^4 <> () or (s^3*j^-1*s^-1)^3 <> () or s^2*j^-1*s^-2*j <> () then return false; fi;

  # check if the generated subgroup acts transitively
  index := Maximum(LargestMovedPoint([s, j]), 1);
  return IsTransitive(Group(s,j), [1..index]);
end);

InstallMethod(CosetActionFromGenerators, [IsRectangularTable], function(gens)
  local A, F2, S, T, SL2Z, gen_words, decomp, coset_table, i;
  if IsEmpty(gens) or not IsMatrix(gens[1]) then
    Error("<gens> needs to be a non-empty list of matrices in SL(2,Z)");
  fi;
  for A in gens do
    if not A in SL(2,Integers) then
      Error("<gens> needs to be a non-empty list of matrices in SL(2,Z)");
    fi;
  od;


  # GAP can only enumerate cosets in the setting of finitely presented groups,
  # so we need to use a presentation of SL(2,Z) and write every generator matrix
  # as a word in S and T
  F2 := FreeGroup("S", "T");
  S := F2.1; T := F2.2;
  SL2Z := F2 / ParseRelators([S,T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");

  # write every generator matrix as a word in S and T
  gen_words := [];
  for i in [1..Length(gens)] do
    decomp := STDecomposition(gens[i]);
    gen_words[i] := ElementOfFpGroup(FamilyObj(SL2Z.1), ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(decomp)));
  od;

  # enumerate the cosets and find a coset table.
  coset_table := List(CosetTable(SL2Z, Subgroup(SL2Z, gen_words)), PermList);

  return [coset_table[1], coset_table[3]];
end);

InstallMethod(STDecomposition, [IsMatrix], function(M)
  local MatS, MatT, F2, S, T, SL2Z, decomposition, k;

  if not M in SL(2,Integers) then
    Error("<M> needs to be in SL(2,Z)");
  fi;

  MatS := [[0,-1],[1,0]]; MatT := [[1,1],[0,1]];

  F2 := FreeGroup("S", "T");
  S := F2.1; T := F2.2;
  SL2Z := F2 / ParseRelators([S,T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");
  S := SL2Z.1; T := SL2Z.2;

  decomposition := One(SL2Z);

  while M[2][1] <> 0 do
    k := QuoInt(M[2][2], M[2][1]);
    decomposition := S^-1 * T^k * decomposition;
    M := M * MatT^-k * MatS;
  od;

  # now M[2][1] = 0 and since det M = 1, we have M = +/- T^r
  # where r = M[1][2]
  if M[1][1] = 1 then
    decomposition := T^(M[1][2]) * decomposition;
  else
    decomposition := S^2 * T^(-M[1][2]) * decomposition;
  fi;

  return decomposition;
end);

InstallMethod(RTDecomposition, [IsMatrix], function(M)
  local w, F2ST, F2RT;

  w := STDecomposition(M);

  F2ST := FreeGroup("S", "T");
  F2RT := FreeGroup("R", "T");

  w := ObjByExtRep(FamilyObj(F2ST.1), ExtRepOfObj(w));
  return MappedWord(w, [F2ST.1, F2ST.2], [F2RT.1*F2RT.2^-1*F2RT.1, F2RT.2]);
end);

InstallMethod(SJDecomposition, [IsMatrix], function(M)
  local w, F2ST, F2SJ;

  w := STDecomposition(M);

  F2ST := FreeGroup("S", "T");
  F2SJ := FreeGroup("S", "J");

  w := ObjByExtRep(FamilyObj(F2ST.1), ExtRepOfObj(w));
  return MappedWord(w, [F2ST.1, F2ST.2], [F2SJ.1, F2SJ.2^-1*F2SJ.1^-1]);
end);

InstallMethod(STDecompositionAsList, [IsMatrix], function(M)
  local MatS, MatT, decomposition, k;

  if not M in SL(2,Integers) then
    Error("<M> needs to be in SL(2,Z)");
  fi;

  MatS := [[0,-1],[1,0]]; MatT := [[1,1],[0,1]];

  decomposition := [];

  while M[2][1] <> 0 do
    k := QuoInt(M[2][2], M[2][1]);
    decomposition := Concatenation([["S", -1], ["T", k]], decomposition);
    M := M * MatT^-k * MatS;
  od;

  # now M[2][1] = 0 and since det M = 1, we have M = +/- T^r
  # where r = M[1][2]
  if M[1][1] = 1 then
    decomposition := Concatenation([["T", M[1][2]]], decomposition);
  else
    decomposition := Concatenation([["S", 2], ["T", -M[1][2]]], decomposition);
  fi;

  return decomposition;
end);

InstallMethod(SAction, [IsModularSubgroup], function(G)
  return G!.s;
end);

InstallMethod(TAction, [IsModularSubgroup], function(G)
  return G!.t;
end);

InstallMethod(RAction, [IsModularSubgroup], function(G)
  return G!.r;
end);

InstallMethod(JAction, [IsModularSubgroup], function(G)
  return G!.j;
end);

InstallMethod(CosetActionOf, [IsMatrix, IsModularSubgroup], function(A, G)
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

InstallMethod(Index, "for a modular subgroup", [IsModularSubgroup], function(G)
  local s, t;
  s := SAction(G);
  t := TAction(G);
  return Maximum(LargestMovedPoint([s, t]), 1);
end);

InstallMethod(IsCongruence, "for a modular subgroup", [IsModularSubgroup], function(G)
  # This is an implementation of an algorithm described by Thomas Hamilton
  # and David Loeffler (https://doi.org/10.1112/S1461157013000338)

  local s, t, r, L, N, e, m, a, b, c, d, p, q, u;

  s := SAction(G);
  t := TAction(G);
  r := s^2*t*s^-1*t;

  L := GeneralizedLevel(G);
  if IsElementOf([[-1,0],[0,-1]], G) then
    N := L;
  else
    N := 2 * L;
  fi;

  e := 1;
  while RemInt(N, 2) <> 1 do
    e := e * 2;
    N := N / 2;
  od;
  m := N;
  N := e*m;

  if e = 1 then # N is odd
    return (r^2*t^(-(1/2 mod N)))^3 = ();
  elif m = 1 then # N is a power of 2
    q := t^20*r^(1/5 mod N)*t^-4*r^-1;
    return (t*r^-1*t)^-1*q*t*r^-1*t = q^-1 and
           q^-1*r*q = r^25 and
           (q*r^5*t*r^-1*t)^3 = (t*r^-1*t)^2;
  else
    c := ChineseRem([e, m], [0, 1]) mod N;
    d := ChineseRem([e, m], [1, 0]) mod N;
    a := t^c; b := r^c;
    p := t^d; q := r^d;
    u := p^20*q^(1/5 mod e)*p^-4*q^-1;
    return a*q*a^-1*q^-1 = () and
           (a*b^-1*a)^4 = () and
           (a*b^-1*a)^2 = (b^-1*a)^3 and
           (a*b^-1*a)^2 = (b^2*a^(-(1/2 mod m)))^3 and
           (p*q^-1*p)^-1*u*(p*q^-1*p) = u^-1 and
           u^-1*q*u = q^25 and
           (p*q^-1*p)^2 = (u*q^5*p*q^-1*p)^3;
  fi;
end);

InstallMethod(RightCosetRepresentatives, [IsModularSubgroup], function(G)
  local s, t, F2, SL2Z, S, T, index, H, coset_table;
  s := SAction(G);
  t := TAction(G);

  F2 := FreeGroup("S", "T");
  S := F2.1;
  T := F2.2;
  SL2Z := F2 / ParseRelators([S, T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");

  index := Index(G);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(SL2Z), coset_table);

  # A slicker way to do this without having to explicitly mention a coset table
  # is the following:
  #hom := GroupHomomorphismByImagesNC(SL2Z, Group([s,t]), [SL2Z.1,SL2Z.2], [s,t]);
  #H := PreImage(hom, Stabilizer(Image(hom), 1));
  # Unfortunately though, computing the stabilizer becomes really slow for large (>1000) index

  return AsList(RightTransversal(SL2Z, H));
end);

# The generalized level of G is defined to be the smallest n for which Deficiency(G, n) is minimal.
# We know that Deficiency(G, 2*l) (where l is the Wohlfahrt level) is minimal
InstallMethod(GeneralizedLevel, [IsModularSubgroup], function(G)
  local l, d1, d2;

  l := WohlfahrtLevel(G);
  d1 := Deficiency(G, l);
  d2 := Deficiency(G, 2*l);

  if d1 <= d2 then
    return l;
  else
    return 2*l;
  fi;
end);

# The Wohlfahrt level is the least common multiple of all cusp widths. Thus it only depends on the image of G in PSL(2,Z).
InstallMethod(WohlfahrtLevel, [IsModularSubgroup], function(G)
  local s, t, i, ind, orbits, plist, p, k;
  s := SAction(G);
  t := TAction(G);
  i := s^2;

  if 1^i = 1 then # -1 in G
    return Order(t);
  fi;

  ind := LargestMovedPoint([i, t]);
  orbits := [];
  for k in [1..ind] do
    Add(orbits, Set(Cycle(i, [1..ind], k)));
  od;
  orbits := Set(orbits);
  plist := [];
  for k in [1..Length(orbits)] do
    Add(plist, Position(orbits, OnSets(orbits[k], t)));
  od;
  p := PermList(plist);
  return Order(p);
end);

# The congruence level is the least smallest n such that G contains the principal congruence subgroup \Gamma(n).
InstallMethod(CongruenceLevel, [IsModularSubgroup], function(G)
  local l, one;
  
  if not IsCongruence(G) then
    Print("This group is not a congruence subgroup.\n");
    return fail;
  fi;

  # if G is congruence, then the congruence level of G equals either l or 2*l, where l is the Wohlfahrt level
  # reference: 'Lifts of projective congruence groups', Kiming, Schütt, Verrill; Journal of the LMS, 2011

  l := WohlfahrtLevel(G);
  one := [[1,0],[0,1]];
  if (-one in G) or Deficiency(G, l) = 1 then
    return l;
  else
    return 2*l;
  fi;
end);

InstallOtherMethod(GeneratorsOfGroup, [IsModularSubgroup], function(G)
  local s, t, F2, S, T, SL2Z, coset_table, H, index;

  s := SAction(G);
  t := TAction(G);

  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  SL2Z := F2 / ParseRelators([S,T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");

  index := Index(G);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(SL2Z), coset_table);

  return GeneratorsOfGroup(H);
end);

InstallMethod(MatrixGeneratorsOfGroup, [IsModularSubgroup], function(G)
  local gens, F2, MatS, MatT;
  gens := ShallowCopy(GeneratorsOfGroup(G));
  F2 := FreeGroup("S", "T");
  MatS := [[0,-1],[1,0]];
  MatT := [[1,1],[0,1]];
  Apply(gens, w -> ObjByExtRep(FamilyObj(F2.1), ExtRepOfObj(w)));
  Apply(gens, w -> MappedWord(w, [F2.1, F2.2], [MatS, MatT]));
  return gens;
end);

InstallMethod(IsElementOf, [IsMatrix, IsModularSubgroup], function(A, G)
  local p;

  if not A in SL(2,Integers) then
    Error("<A> has to be an element of SL(2,Z)");
  fi;

  # A matrix A is in the group G if and only if it does not permute the identity coset.
  p := CosetActionOf(A, G);
  return 1^p = 1;
end);
InstallMethod(\in, "for a finite-index subgroup of SL(2,Z)", [IsMatrix, IsModularSubgroup], 10000, function(A, G)
  return IsElementOf(A, G);
end);

InstallMethod(IsSubset, "for two finite-index subgroups of SL(2,Z)", [IsModularSubgroup, IsModularSubgroup], function(H, G)
  local gens, result, g;
  gens := MatrixGeneratorsOfGroup(H);
  result := true;
  for g in gens do
    result := result and (g in G);
  od;
  return result;
end);

InstallMethod(\=, "for two finite-index subgroups of SL(2,Z)", [IsModularSubgroup, IsModularSubgroup], function(G, H)
  return (Index(G) = Index(H)) and IsSubgroup(G, H);
end);

InstallMethod(CuspWidth, [IsRat, IsModularSubgroup], function(c, G)
  return CuspWidth(c, Projection(G));
end);
InstallOtherMethod(CuspWidth, [IsInfinity, IsModularSubgroup], function(c, G)
  return CuspWidth(c, Projection(G));
end);

InstallMethod(CuspsEquivalent, [IsRat, IsRat, IsModularSubgroup], function(c1, c2, G)
  return CuspsEquivalent(c1, c2, Projection(G));
end);
InstallOtherMethod(CuspsEquivalent, [IsRat, IsInfinity, IsModularSubgroup], function(c1, c2, G)
  return CuspsEquivalent(c1, c2, Projection(G));
end);
InstallOtherMethod(CuspsEquivalent, [IsInfinity, IsRat, IsModularSubgroup], function(c1, c2, G)
  return CuspsEquivalent(c2, c1, G);
end);
InstallOtherMethod(CuspsEquivalent, [IsInfinity, IsInfinity, IsModularSubgroup], function(c1, c2, G)
  return true;
end);

InstallMethod(Cusps, [IsModularSubgroup], function(G)
  return Cusps(Projection(G));
end);

InstallMethod(CosetRepresentativeOfCusp, [IsRat, IsModularSubgroup], function(c, G)
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

InstallMethod(MoebiusTransformation, [IsMatrix, IsRat], function(A, r)
  local a, b, c, d, p, q;
  if not A in SL(2,Integers) then
    Error("<A> needs to be a matrix in SL(2,Z).");
  fi;

  a := A[1][1]; b := A[1][2];
  c := A[2][1]; d := A[2][2];
  p := NumeratorRat(r);
  q := DenominatorRat(r);
  if p = -d and q = c then
    return infinity;
  fi;

  return (a*r+b)/(c*r+d);
end);
InstallOtherMethod(MoebiusTransformation, [IsMatrix, IsInfinity], function(A, i)
  if not A in SL(2,Integers) then
    Error("<A> needs to be a matrix in SL(2,Z).");
  fi;
  if A[2][1] = 0 then
    return infinity;
  fi;
  return A[1][1]/A[2][1];
end);

InstallMethod(IndexModN, [IsModularSubgroup, IsPosInt], function(G, N)
  local gens, SL2Zn, H;
  if N = 1 then return 1; fi;
  gens := ShallowCopy(MatrixGeneratorsOfGroup(G));
  Apply(gens, M ->
    [[ZmodnZObj(M[1][1], N), ZmodnZObj(M[1][2], N)],
     [ZmodnZObj(M[2][1], N), ZmodnZObj(M[2][2], N)]]
  );
  SL2Zn := SL(2, Integers mod N);
  H := Subgroup(SL2Zn, gens);
  return Index(SL2Zn, H);
end);

InstallMethod(Deficiency, [IsModularSubgroup], function(G)
  if IsCongruence(G) then return 1; fi;
  return Index(G) / IndexModN(G, 2 * GeneralizedLevel(G));
end);

InstallMethod(Deficiency, [IsModularSubgroup, IsPosInt], function(G, N)
  if HasGeneralizedLevel(G) then
    # if G is congruence, then the congruence level of G equals either l or 2*l, where l is the Wohlfahrt level
    # reference: 'Lifts of projective congruence groups', Kiming, Schütt, Verrill; Journal of the LMS, 2011
    if RemInt(N, 2*GeneralizedLevel(G)) = 0 and IsCongruence(G) then return 1; fi;
  fi;
  return Index(G) / IndexModN(G, N);
end);

InstallMethod(Projection, [IsModularSubgroup], function(G)
  local s, t, i, ind, orbits, plist, p, k, qlist, q;
  s := SAction(G);
  t := TAction(G);
  i := s^2;

  if 1^i = 1 then # -1 in G
    return ProjectiveModularSubgroup(s,t);
  fi;

  # calculate the action of T on the cosets in PSL(2,Z)
  ind := LargestMovedPoint([i, t]);
  orbits := [];
  # group those cosets together which are identified by -1
  for k in [1..ind] do
    Add(orbits, Set(Cycle(i, [1..ind], k)));
  od;
  orbits := Set(orbits);
  plist := [];
  for k in [1..Length(orbits)] do
    Add(plist, Position(orbits, OnSets(orbits[k], t)));
  od;
  p := PermList(plist);


  # calculate the action of S on the cosets in PSL(2,Z)
  ind := LargestMovedPoint([i, s]); # is this necessary??
  orbits := [];
  for k in [1..ind] do
    Add(orbits, Set(Cycle(i, [1..ind], k)));
  od;
  orbits := Set(orbits);
  qlist := [];
  for k in [1..Length(orbits)] do
    Add(qlist, Position(orbits, OnSets(orbits[k], s)));
  od;
  q := PermList(qlist);


  return ProjectiveModularSubgroup(q, p);
end);

InstallMethod(Conjugate, [IsModularSubgroup, IsMatrix], function(G, A)
  local a, s, t;
  a := CosetActionOf(A, G);
  s := SAction(G);
  t := TAction(G);
  return ModularSubgroup(a^-1*s*a, a^-1*t*a); # this is A^-1*G*A
  # note that if G = (\sigma_S, \sigma_T) then A^-1*G*A = (\sigma_A*\sigma_S*\sigma_A^-1, \sigma_A*\sigma_T*\sigma_A^-1)
  # but GAP multiplies permutations the other way around
end);

InstallMethod(NormalCore, [IsModularSubgroup], function(G)
  local s, t, F2, S, T, SL2Z, index, m, Sd, core, coset_table;

  s := SAction(G);
  t := TAction(G);

  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  SL2Z := F2 / ParseRelators([S,T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");
  S := SL2Z.1;
  T := SL2Z.2;

  index := Index(G);
  Sd := SymmetricGroup(index);

  m := GroupHomomorphismByImagesNC(SL2Z, Sd, [S,T], [s,t]);
  core := Kernel(m);
  coset_table := CosetTableBySubgroup(SL2Z, core);
  return ModularSubgroup(PermList(coset_table[1]), PermList(coset_table[3]));
end);

InstallMethod(QuotientByNormalCore, [IsModularSubgroup], function(G)
  local NC, s, t, F2, S, T, SL2Z, index, coset_table, H;

  NC := NormalCore(G);
  s := SAction(NC);
  t := TAction(NC);

  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  SL2Z := F2 / ParseRelators([S,T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");

  index := Index(NC);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(SL2Z), coset_table);

  return FactorGroupNC(SL2Z, H);
end);

InstallMethod(AssociatedCharacterTable, [IsModularSubgroup], function(G)
  local F;
  F := QuotientByNormalCore(G);
  return CharacterTable(F);
end);

InstallMethod(Genus, [IsModularSubgroup], function(G)
  return Genus(Projection(G));
end);

InstallMethod(IO_Pickle, "for a modular subgroup", [IsFile, IsModularSubgroup], function(f, G)
  local nr;
  nr := IO_AddToPickled(G); # this should always be false (?)
  if IO_Write(f, "MODS") = fail then IO_FinalizePickled(); return IO_Error; fi;
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

  if HasMatrixGeneratorsOfGroup(G) then
    if IO_Pickle(f, MatrixGeneratorsOfGroup(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

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

  if HasDeficiency(G) then
    if IO_Pickle(f, Deficiency(G)) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  else
    if IO_Pickle(f, fail) = IO_Error then IO_FinalizePickled(); return IO_Error; fi;
  fi;

  IO_FinalizePickled();
  return IO_OK;
end);

IO_Unpicklers.MODS := function(f)
  local H, s, t, cong, cusps, level, matrixgen, normal_core, quotient_nc, ctbl, genus, deficiency;
  s := IO_Unpickle(f);
  if s = IO_Error then return IO_Error; fi;
  t := IO_Unpickle(f);
  if t = IO_Error then return IO_Error; fi;
  H := ModularSubgroup(s, t);
  IO_AddToUnpickled(H);

  cong := IO_Unpickle(f);
  if cong = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  cusps := IO_Unpickle(f);
  if cusps = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  level := IO_Unpickle(f);
  if level = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  matrixgen := IO_Unpickle(f);
  if matrixgen = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  normal_core := IO_Unpickle(f);
  if normal_core = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  quotient_nc := IO_Unpickle(f);
  if quotient_nc = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  ctbl := IO_Unpickle(f);
  if ctbl = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  genus := IO_Unpickle(f);
  if genus = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;
  deficiency := IO_Unpickle(f);
  if deficiency = IO_Error then IO_FinalizeUnpickled(); return IO_Error; fi;

  if cong <> fail then SetIsCongruence(H, cong); fi;
  if cusps <> fail then SetCusps(H, cusps); fi;
  if level <> fail then SetGeneralizedLevel(H, level); fi;
  if matrixgen <> fail then SetMatrixGeneratorsOfGroup(H, matrixgen); fi;
  if normal_core <> fail then SetNormalCore(H, normal_core); fi;
  if quotient_nc <> fail then SetQuotientByNormalCore(H, quotient_nc); fi;
  if ctbl <> fail then SetAssociatedCharacterTable(H, ctbl); fi;
  if genus <> fail then SetGenus(H, genus); fi;
  if deficiency <> fail then SetDeficiency(H, deficiency); fi;

  IO_FinalizeUnpickled();
  return H;
end;

InstallMethod(PrintObj, "for modular subgroups", [IsModularSubgroup], function(G)
  Print("ModularSubgroup(",
          "\n  S : ", SAction(G),
          "\n  T : ", TAction(G),
          "\n  R : ", RAction(G),
          "\n  J : ", JAction(G),
        " )");
end);

InstallMethod(ViewObj, "for modular subgroups", [IsModularSubgroup], function(G)
  Print("<modular subgroup of index ", Index(G),">");
end);
