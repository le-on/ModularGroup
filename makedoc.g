LoadPackage( "AutoDoc" );
AutoDoc( rec(
  autodoc := true,
  scaffold := true,
  gapdoc := rec(
    scan_dirs := ["lib"],
    files := ["lib/ModularSubgroups.gd", "lib/ModularSubgroups.gi"]
  )
) );
QUIT;
