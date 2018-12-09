DeclareCategory("IsModularSubgroup", IsMatrixGroup);
DeclareRepresentation("IsModularSubgroupRepresentation", IsComponentObjectRep, ["s", "t", "r", "j"]);
DeclareSynonym("IsDefaultModularSubgroup", IsModularSubgroup and IsModularSubgroupRepresentation);


DeclareAttribute("Index", IsModularSubgroup);
DeclareAttribute("IsCongruence", IsModularSubgroup);
DeclareAttribute("Cusps", IsModularSubgroup);
DeclareAttribute("RightCosetRepresentatives", IsModularSubgroup);
DeclareAttribute("GeneralizedLevel", IsModularSubgroup);
DeclareAttribute("GeneratorsOfGroup", IsModularSubgroup);
DeclareAttribute("MatrixGeneratorsOfGroup", IsModularSubgroup);
DeclareAttribute("NormalCore", IsModularSubgroup);
DeclareAttribute("QuotientByNormalCore", IsModularSubgroup);
DeclareAttribute("AssociatedCharacterTable", IsModularSubgroup);
DeclareAttribute("Genus", IsModularSubgroup);

DeclareOperation("DefinesCosetActionST", [IsPerm, IsPerm]);
DeclareOperation("DefinesCosetActionRT", [IsPerm, IsPerm]);
DeclareOperation("DefinesCosetActionSJ", [IsPerm, IsPerm]);
DeclareOperation("ModularSubgroup", [IsPerm, IsPerm]);
DeclareOperation("ModularSubgroupST", [IsPerm, IsPerm]);
DeclareOperation("ModularSubgroupRT", [IsPerm, IsPerm]);
DeclareOperation("ModularSubgroupSJ", [IsPerm, IsPerm]);
DeclareOperation("SAction", [IsModularSubgroup]);
DeclareOperation("TAction", [IsModularSubgroup]);
DeclareOperation("RAction", [IsModularSubgroup]);
DeclareOperation("JAction", [IsModularSubgroup]);
DeclareOperation("CosetActionOf", [IsMatrix, IsModularSubgroup]);
DeclareOperation("CosetActionFromGenerators", [IsRectangularTable]);
DeclareOperation("STDecomposition", [IsMatrix]);
DeclareOperation("RTDecomposition", [IsMatrix]);
DeclareOperation("SJDecomposition", [IsMatrix]);
DeclareOperation("STDecompositionAsList", [IsMatrix]);
DeclareOperation("IsElementOf", [IsMatrix, IsModularSubgroup]);
DeclareOperation("CuspWidth", [IsRat, IsModularSubgroup]);
DeclareOperation("CuspsEquivalent", [IsRat, IsRat, IsModularSubgroup]);
DeclareOperation("IndexModN", [IsModularSubgroup, IsPosInt]);
DeclareOperation("Deficiency", [IsModularSubgroup, IsPosInt]);
DeclareOperation("Projection", [IsModularSubgroup]);

DeclareOperation("MoebiusTransformation", [IsMatrix, IsRat]);
