InstallMethod(ModularSubgroup, [IsPerm, IsPerm], function(sp, tp)
  local G, type, tab, index, l1, l2, l3, l4;

  if not DefinesCosetAction(sp, tp) then
    Error("<s> and <t> do not describe the action of the generators S and T on the cosets of a finite-index subgroup of SL(2,Z)");
  fi;

  type := NewType(FamilyObj(One(SL(2,Integers))),
    IsObject and
    #IsMatrixGroup and
    IsAttributeStoringRep and
    IsComponentObjectRep and
    IsFinitelyGeneratedGroup and # finite-index subgroups of finitely generated groups are finitely generated
    IsDefaultModularSubgroup);

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

InstallMethod(DefinesCosetAction, [IsPerm, IsPerm], function(s, t)
  local index;

  # check relations
  if s^4 <> () or (s^3*t)^3 <> () or s^2*t*s^-2*t^-1 <> () then return false; fi;

  # check if the generated subgroup acts transitively
  # this check can be quite costly, so we do it last
  index := Maximum(LargestMovedPoint([s, t]), 1);
  return IsTransitive(Group(s,t), [1..index]);
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

InstallMethod(SAction, [IsModularSubgroup], function(G)
  return G!.s;
end);

InstallMethod(TAction, [IsModularSubgroup], function(G)
  return G!.t;
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

InstallMethod(GeneralizedLevel, [IsModularSubgroup], function(G)
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

InstallMethod(GeneratorsOfGroup, [IsModularSubgroup], function(G)
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
  gens := ShallowCopy(MatrixGeneratorsOfGroup(G));
  Apply(gens, M ->
    [[ZmodnZObj(M[1][1], N), ZmodnZObj(M[1][2], N)],
     [ZmodnZObj(M[2][1], N), ZmodnZObj(M[2][2], N)]]
  );
  SL2Zn := SL(2, Integers mod N);
  H := Subgroup(SL2Zn, gens);
  return Index(SL2Zn, H);
end);

InstallMethod(Deficiency, [IsModularSubgroup, IsPosInt], function(G, N)
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

InstallMethod(PrintObj, "for modular subgroups", [IsModularSubgroup], function(G)
  Print("ModularSubgroup( ", SAction(G), ", ", TAction(G)," )");
end);

InstallMethod(ViewObj, "for modular subgroups", [IsModularSubgroup], function(G)
  Print("<modular subgroup of index ", Index(G),">");
end);
