DeclareCategory("IsProjectiveModularSubgroup", IsFinitelyGeneratedGroup);
DeclareRepresentation("IsProjectiveModularSubgroupRepresentation", IsComponentObjectRep, ["s", "t"]);
DeclareSynonym("IsDefaultProjectiveModularSubgroup", IsProjectiveModularSubgroup and IsProjectiveModularSubgroupRepresentation);


DeclareAttribute("Index", IsProjectiveModularSubgroup);
DeclareAttribute("IsCongruence", IsProjectiveModularSubgroup);
DeclareAttribute("Cusps", IsProjectiveModularSubgroup);
DeclareAttribute("RightCosetRepresentatives", IsProjectiveModularSubgroup);
DeclareAttribute("GeneralizedLevel", IsProjectiveModularSubgroup);
DeclareAttribute("GeneratorsOfGroup", IsProjectiveModularSubgroup);
DeclareAttribute("NormalCore", IsProjectiveModularSubgroup);
DeclareAttribute("QuotientByNormalCore", IsProjectiveModularSubgroup);
DeclareAttribute("AssociatedCharacterTable", IsProjectiveModularSubgroup);
DeclareAttribute("Genus", IsProjectiveModularSubgroup);

DeclareOperation("DefinesProjectiveCosetAction", [IsPerm, IsPerm]);
DeclareOperation("ProjectiveModularSubgroup", [IsPerm, IsPerm]);
DeclareOperation("SAction", [IsProjectiveModularSubgroup]);
DeclareOperation("TAction", [IsProjectiveModularSubgroup]);
DeclareOperation("CosetActionOf", [IsMatrix, IsProjectiveModularSubgroup]);
DeclareOperation("IsElementOf", [IsMatrix, IsProjectiveModularSubgroup]);
DeclareOperation("CuspWidth", [IsRat, IsProjectiveModularSubgroup]);
DeclareOperation("CuspsEquivalent", [IsRat, IsRat, IsProjectiveModularSubgroup]);
DeclareOperation("CosetRepresentativeOfCusp", [IsRat, IsProjectiveModularSubgroup]);
DeclareOperation("LiftToSL2ZEven", [IsProjectiveModularSubgroup]);
DeclareOperation("LiftToSL2ZOdd", [IsProjectiveModularSubgroup]);
DeclareOperation("IndexModN", [IsProjectiveModularSubgroup, IsPosInt]);
DeclareOperation("Deficiency", [IsProjectiveModularSubgroup, IsPosInt]);
