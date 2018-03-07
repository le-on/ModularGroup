DeclareCategory("IsModularSubgroup", IsObject);
DeclareRepresentation("IsModularSubgroupRepresentation", IsComponentObjectRep, ["s", "t"]);
DeclareSynonym("IsDefaultModularSubgroup", IsModularSubgroup and IsModularSubgroupRepresentation);


DeclareAttribute("Index", IsModularSubgroup);
DeclareAttribute("IsCongruenceSubgroup", IsModularSubgroup);
DeclareAttribute("Cusps", IsModularSubgroup);
DeclareAttribute("RightCosetRepresentatives", IsModularSubgroup);
DeclareAttribute("GeneralizedLevel", IsModularSubgroup);
DeclareAttribute("GeneratorsOfGroup", IsModularSubgroup);

DeclareOperation("DefinesCosetAction", [IsPerm, IsPerm]);
DeclareOperation("ModularSubgroup", [IsPerm, IsPerm]);
DeclareOperation("SAction", [IsModularSubgroup]);
DeclareOperation("TAction", [IsModularSubgroup]);
DeclareOperation("CosetActionOf", [IsMatrix, IsModularSubgroup]);
DeclareOperation("CosetActionFromGenerators", [IsRectangularTable]);
DeclareOperation("STDecomposition", [IsMatrix]);
DeclareOperation("IsElementOf", [IsMatrix, IsModularSubgroup]);
DeclareOperation("CuspWidth", [IsRat, IsModularSubgroup]);
DeclareOperation("CuspsEquivalent", [IsRat, IsRat, IsModularSubgroup]);
DeclareOperation("IndexModN", [IsModularSubgroup, IsPosInt]);
DeclareOperation("Deficiency", [IsModularSubgroup, IsPosInt]);
DeclareOperation("Projection", [IsModularSubgroup]);

DeclareOperation("MoebiusTransformation", [IsMatrix, IsRat]);
