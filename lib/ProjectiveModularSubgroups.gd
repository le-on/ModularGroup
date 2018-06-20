DeclareCategory("IsProjectiveModularSubgroup", IsObject);
DeclareRepresentation("IsProjectiveModularSubgroupRepresentation", IsComponentObjectRep, ["s", "t"]);
DeclareSynonym("IsDefaultProjectiveModularSubgroup", IsProjectiveModularSubgroup and IsProjectiveModularSubgroupRepresentation);


DeclareAttribute("Index", IsProjectiveModularSubgroup);
DeclareAttribute("IsCongruenceSubgroup", IsProjectiveModularSubgroup);
DeclareAttribute("Cusps", IsProjectiveModularSubgroup);
DeclareAttribute("RightCosetRepresentatives", IsProjectiveModularSubgroup);
DeclareAttribute("GeneralizedLevel", IsProjectiveModularSubgroup);
DeclareAttribute("GeneratorsOfGroup", IsProjectiveModularSubgroup);
DeclareAttribute("NormalCore", IsProjectiveModularSubgroup);
DeclareAttribute("QuotientByNormalCore", IsProjectiveModularSubgroup);
DeclareAttribute("AssociatedCharacterTable", IsProjectiveModularSubgroup);

DeclareOperation("DefinesProjectiveCosetAction", [IsPerm, IsPerm]);
DeclareOperation("ProjectiveModularSubgroup", [IsPerm, IsPerm]);
DeclareOperation("SAction", [IsProjectiveModularSubgroup]);
DeclareOperation("TAction", [IsProjectiveModularSubgroup]);
#DeclareOperation("CosetActionFromGenerators", [IsRectangularTable]);
#DeclareOperation("IsElementOf", [IsMatrix, IsProjectiveModularSubgroup]);
DeclareOperation("CuspWidth", [IsRat, IsProjectiveModularSubgroup]);
DeclareOperation("CuspsEquivalent", [IsRat, IsRat, IsProjectiveModularSubgroup]);
DeclareOperation("LiftToSL2ZEven", [IsProjectiveModularSubgroup]);
DeclareOperation("LiftToSL2ZOdd", [IsProjectiveModularSubgroup]);
DeclareOperation("IndexModN", [IsProjectiveModularSubgroup, IsPosInt]);
DeclareOperation("Deficiency", [IsProjectiveModularSubgroup, IsPosInt]);
