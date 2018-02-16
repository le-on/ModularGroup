#! @Chapter The Modular Group and its subgroups
#! This package contains methods for computing with finite-index subgroups of the
#! modular group $\mathrm{SL}(2, \mathbb{Z})$ which are given by a coset permutation
#! representation with respect to the generators
#! @BeginLatexOnly
#! $$ S = \left( {\begin{array}{cc} 0 & -1 \\ 1 & 0 \\ \end{array} } \right)  \quad T = \left( {\begin{array}{cc} 1 & 1 \\ 0 & 1 \\ \end{array} } \right) $$
#! @EndLatexOnly
#! We will call these subgroups 'modular subgroups'.


#! @Section Construction of modular subgroups
#! In this section we describe how to construct modular subgroups from a given
#! coset permutation representation or from a list of generator matrices and
#! some related methods.

#! @Arguments s, t
#! @Returns an object representing a modular subgroup
#! @Label for two permutations
#! @Description
#!  This method constructs a modular subgroup from two given permutations
#!  (provided they describe a coset action).
InstallMethod(ModularSubgroup, [IsPerm, IsPerm], function(sp, tp)
  local G, type, H, tab;

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

  H := Group(sp, tp);
  tab := CosetTableBySubgroup(H, Stabilizer(H, 1));
  StandardizeTable(tab);

  G := Objectify(type, rec(
    s := PermList(tab[1]),
    t := PermList(tab[3])
  ));
  return G;
end);

#! @Arguments gens
#! @Returns an object representing a modular subgroup
#! @Label for a list of matrices
#! @Description
#!  This is another constructor for a modular subgroup, with the difference that it
#!  takes a list of generator matrices as input and calculates the coset graph of
#!  the generated group. One has to be careful though when using this method,
#!  because no check is performed if the generated group has finite index!
#!  Internally, when trying to calculate the coset graph, we just enumerate all
#!  cosets until we are done or some limit fixed is reached. This also exposes
#!  another weakness of this method: The coset enumeration might be very
#!  time-consuming, so constructing modular subgroups from a list of generators
#!  is not always feasible. On the other hand, the clear advantage of
#!  constructing a modular subgroup in this way is that it will alreasy know its
#!  generators. So future computations with this group involving generators will
#!  most likely be faster.
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
  SetGeneratorsOfGroup(G, gens);
  return G;
end);

#! @Arguments s, t
#! @Returns true or false
#! @Label for two permutations
#! @Description
#!  This is an auxiliary method that takes two permutations as input and checks
#!  if they describe the action of the generators $S$ and $T$ on the cosets of
#!  some group. This check is for example performed in the constructor for
#!  a modular subgroup.
InstallMethod(DefinesCosetAction, [IsPerm, IsPerm], function(s, t)
  local index;

  # check relations
  if s^4 <> () or (s^3*t)^3 <> () or s^2*t*s^-2*t^-1 <> () then return false; fi;

  # check if the generated subgroup acts transitively
  # this check can be quite costly, so we do it last
  index := Maximum(LargestMovedPoint([s, t]), 1);
  return IsTransitive(Group(s,t), [1..index]);
end);

#! @Arguments gens
#! @Returns two permutations
#! @Label for a list of generator matrices
#! @Description
#!  Takes a list of generator matrices and calculates the coset permutation
#!  representation of the generated subgroup. The same warning as above applies:
#!  No check is performed if the generated subgroup actually has finite index.
InstallMethod(CosetActionFromGenerators, [IsRectangularTable], function(gens)
  local F2, S, T, SL2Z, gen_words, decomp, coset_table, i;
  if not (IsMatrix(gens[1]) and gens[1] in SL(2,Integers)) then
    Error("<gens> needs to be a non-empty list of matrices in SL(2,Z)");
  fi;


  # GAP can only enumerate cosets in the setting of finitely presented groups, so we must use a presentation
  # of SL(2,Z) and write every generator matrix as a word in S and T.
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

#! @Arguments M
#! @Returns a word in $S$ and $T$
#! @Label for a matrix in SL(2,Z)
#! @Description
#!  Takes a matrix in $\mathrm{SL}(2, \mathbb{Z})$ and decomposes it in a word
#!  in the generators $S$ and $T$.
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

#! @Arguments G
#! @Returns a permutation
#! @Label for a modular subgroup
#! @Description
#!  Takes a modular subgroup and returns a permutation corresponding to the
#!  action of the generator matrix $S$.
InstallMethod(SAction, [IsModularSubgroup], function(G)
  return G!.s;
end);

#! @Arguments G
#! @Returns a permutation
#! @Label for a modular subgroup
#! @Description
#!  Takes a modular subgroup and returns a permutation corresponding to the
#!  action of the generator matrix $T$.
InstallMethod(TAction, [IsModularSubgroup], function(G)
  return G!.t;
end);


#! @Section Computations with modular subgroups
#! In this section we describe the implemented method for computing with
#! modular subgroups.

#! @Arguments G
#! @Returns a natural number
#! @Label for a modular subgroup
#! @Description
#!  Takes a modular subgroup and returns its index in $\mathrm{SL}(2, \mathbb{Z})$.
InstallMethod(Index, "for a modular subgroup", [IsModularSubgroup], function(G)
  local s, t;
  s := SAction(G);
  t := TAction(G);
  return Maximum(LargestMovedPoint([s, t]), 1);
end);

#! @Arguments G
#! @Returns true or false
#! @Label for a modular subgroup
#! @Description
#!  Tests whether a given modular subgroup is a congruence subgroup.
InstallMethod(IsCongruenceSubgroup, "for a modular subgroup", [IsModularSubgroup], function(G)
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

#! @Arguments G
#! @Returns a list of matrices
#! @Label for a modular subgroup
#! @Description
#!  Calculates a list of representatives of the right cosets of a given
#!  modular subgroup.
InstallMethod(RightCosetRepresentatives, [IsModularSubgroup], function(G)
  local s, t, F2, SL2Z, S, T, hom, H;
  s := SAction(G);
  t := TAction(G);

  F2 := FreeGroup("S", "T");
  S := F2.1;
  T := F2.2;
  SL2Z := F2 / ParseRelators([S, T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");

  hom := GroupHomomorphismByImagesNC(SL2Z, Group([s,t]), [SL2Z.1,SL2Z.2], [s,t]);
  H := PreImage(hom, Stabilizer(Image(hom), 1));

  return AsList(RightTransversal(SL2Z, H));
end);

#! @Arguments G
#! @Returns a natural number
#! @Label for a modular subgroup
#! @Description
#!  Computes the generalized level (i.e. the lowest common multiple of all cusp
#!  widths) of a given modular subgroup.
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

#! @Arguments G
#! @Returns
#! @Label for a modular subgroup
#! @Description
#!  Calculates a list of generators for a given modular subgroup.
#!  Note: The returned list might contain redundant generators (or even
#!  duplicates). This calculation involves enumerating the cosets of the given
#!  group and might become very slow for large index.
InstallMethod(GeneratorsOfGroup, [IsModularSubgroup], function(G)
  local s, t, F2, S, T, SL2Z_fp, MatS, MatT, SL2Z, iso, coset_table, H, gens, index;

  s := SAction(G);
  t := TAction(G);
  # Since GAP can only reconstruct generators of a group which is given by a coset graph if the
  # group is a subgroup of a finitely presented group, we need a presentation of SL(2,Z).
  # We will use the following one:
  F2 := FreeGroup("S", "T");
  S := F2.S;
  T := F2.T;
  SL2Z_fp := F2 / ParseRelators([S,T], "S^4, (S^3*T)^3, S^2*T*S^-2*T^-1");
  # "S^4, (S^-1*T^-1)^3*S^-2"

  # SL(2,Z) as matrix group.
  MatS := [[0,-1],[1,0]];
  MatT := [[1,1],[0,1]];
  SL2Z := Group([MatS, MatT]);

  # This is an explicit isomorphism between the above presentation of SL(2,Z) and SL(2,Z) as a matrix group.
  iso := GroupHomomorphismByImagesNC(SL2Z_fp, SL2Z, GeneratorsOfGroup(SL2Z_fp), [MatS,MatT]);

  # In order to compute generators of a subgroup, GAP needs a coset table in terms of generators of the whole
  # group in a specific form. For details, see
  # https://www.gap-system.org/Manuals/doc/ref/chap47.html#X857F239583AFE0B7 and
  # https://www.gap-system.org/Manuals/doc/ref/chap47.html#X7F7F31E47D7F6EF8
  index := Index(G);
  coset_table := [ListPerm(s, index), ListPerm(s^-1, index), ListPerm(t, index), ListPerm(t^-1, index)];
  H := SubgroupOfWholeGroupByCosetTable(FamilyObj(SL2Z_fp), coset_table);

  # The method 'Apply' used below requires its argument to be a mutable list, but 'GeneratorsOfGroup' returns
  # an immutable list, so we need to make a mutable copy of it.
  gens := ShallowCopy(GeneratorsOfGroup(H));

  # Now we apply the isomorphism defined above to get the generators as matrices in SL(2,Z) rather than
  # words in S and T.
  Apply(gens, x -> Image(iso, x));

  return gens;
end);

#! @Arguments A, G
#! @Returns true or false
#! @Label for a matrix in SL(2,Z) and a modular subgroup
#! @Description
#!  This is a membership test for modular subgroups given by a coset permutation
#!  representation.
InstallMethod(IsElementOf, [IsMatrix, IsModularSubgroup], function(A, G)
  local s, t, decomp, rep, p, i;

  if not A in SL(2,Integers) then
    Error("<A> has to be an element of SL(2,Z)");
  fi;

  s := SAction(G);
  t := TAction(G);

  # A matrix A is in the group if and only if it does not permute the identity coset.
  # Hence, we write A as a word in S and T (the action of which we know)
  # and compute from that the action of A.

  decomp := STDecomposition(A);
  rep := ExtRepOfObj(decomp);

  if Length(rep) = 0 then
    return s = () and t = ();
  fi;

  p := ();
  for i in [1, 3 .. Length(rep)-1] do
    if rep[i] = 1 then
      p := p * s^rep[i+1]; # the action is from the right
    else
      p := p * t^rep[i+1];
    fi;
  od;
  return 1^p = 1;
end);

#! @Arguments c, G
#! @Returns a natural number
#! @Label for a rational number or infinity and a modular subgroup
#! @Description
#!  Calculates the width of $c$ with respect to a given modular subgroup, i.e.
#!  the smallest $k$ such that $\pm gT^{k}g^{-1} \in G$ where $g \in \mathrm{SL}(2, \mathbb{Z})$
#!  such that $g\infty = c$.
InstallMethod(CuspWidth, [IsRat, IsModularSubgroup], function(c, G)
  return CuspWidth(c, Projection(G));
end);
InstallOtherMethod(CuspWidth, [IsInfinity, IsModularSubgroup], function(c, G)
  return CuspWidth(c, Projection(G));
end);


#! @Arguments c1, c2, G
#! @Returns true or false
#! @Label for two cusps (i.e. rational numbers or infinity) and a modular subgroup
#! @Description
#!  Checks if two cusps are equivalent with respect to a given modular subgroup.
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

#! @Arguments G
#! @Returns a list of cusps
#! @Label for a modular subgroup
#! @Description
#!  Calculates a list of inequivalent cusp representative for a given modular subgroup.
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

#! @Arguments G, N
#! @Returns a natural number
#! @Label for a modular subgroup and a natural number > 1
#! @Description
#!  This method computes the index of the image of $G$ in $\mathrm{SL}(2, \mathbb{Z}/N\mathbb{Z})$
#!  under the projection
#!  $$\pi_N : \mathrm{SL}(2,\mathbb{Z}) \rightarrow \mathrm{SL}(2, \mathbb{Z}/N\mathbb{Z})$$
InstallMethod(IndexModN, [IsModularSubgroup, IsPosInt], function(G, N)
  local gens, SL2Zn, H;
  gens := ShallowCopy(GeneratorsOfGroup(G));
  Apply(gens, M ->
    [[ZmodnZObj(M[1][1], N), ZmodnZObj(M[1][2], N)],
     [ZmodnZObj(M[2][1], N), ZmodnZObj(M[2][2], N)]]
  );
  SL2Zn := SL(2, Integers mod N);
  H := Subgroup(SL2Zn, gens);
  return Index(SL2Zn, H);
end);

#! @Arguments G, N
#! @Returns a natural number
#! @Label for a modular subgroup and a natural number > 1
#! @Description
#!  This method calculates the so-called deficiency $f_N$ of a modular subgroup,
#!  i.e. the index $[ \Gamma(N) : \Gamma(N) \cap G ]$.
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

InstallMethod(PrintObj, "for modular subgroups", [IsModularSubgroup], function(G)
  Print("ModularSubgroup( ", SAction(G), ", ", TAction(G)," )");
end);

InstallMethod(ViewObj, "for modular subgroups", [IsModularSubgroup], function(G)
  Print("<modular subgroup of index ", Index(G),">");
end);
