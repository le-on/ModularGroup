ReadPackage("ModularGroup/lib/ModularSubgroups.gi");
ReadPackage("ModularGroup/lib/ProjectiveModularSubgroups.gi");
if TestPackageAvailability("Congruence", ">=1.1.1") <> fail then
  ReadPackage("ModularGroup/lib/misc.gi");
fi;
