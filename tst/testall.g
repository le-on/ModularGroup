LoadPackage("ModularGroup");

TestDirectory(DirectoriesPackageLibrary( "ModularGroup", "tst" ),
  rec(exitGAP     := true,
      testOptions := rec(compareFunction := "uptowhitespace") ) );

FORCE_QUIT_GAP(1);
