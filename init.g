ReadPackage("ModularGroup/lib/ModularSubgroups.gd");
ReadPackage("ModularGroup/lib/ProjectiveModularSubgroups.gd");
if TestPackageAvailability("Congruence", ">=1.1.1") <> fail then
  ReadPackage("ModularGroup/lib/misc.gd");
fi;
