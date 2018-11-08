InstallMethod(LiftToSL2Z, [IsMatrix, IsPosInt], function(a, b)
  Info(InfoWarning, 1, "The 'Congruence' package is not installed. You need to install it to use this function.");
end);
InstallOtherMethod(LiftToSL2Z, [IsMatrixGroup, IsPosInt], function(a, b)
  Info(InfoWarning, 1, "The 'Congruence' package is not installed. You need to install it to use this function.");
end);
