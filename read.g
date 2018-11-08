ReadPackage("ModularGroup/lib/ModularSubgroups.gi");
ReadPackage("ModularGroup/lib/ProjectiveModularSubgroups.gi");
if TestPackageAvailability("Congruence", ">=1.1.1") = fail then
  Info(InfoWarning, 1, "The function 'LiftToSL2Z' relies on the 'Congruence' package which seems to be missing from your GAP installation. If you want to use this function, you need to install the 'Congruence' package.");
  ReadPackage("ModularGroup/lib/misc_no_cong.gi");
else
  ReadPackage("ModularGroup/lib/misc.gi");
fi;
