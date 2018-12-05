################################################################################
##                                                                            ##
##  PackageInfo.g        The `ModularGroup' package                 Luca Junk ##
##                                                                            ##
################################################################################

SetPackageInfo( rec(

PackageName := "ModularGroup",
Subtitle := "Finite-index subgroups of (P)SL(2,Integers)",
Version := "1.0.0",
Date := "05/12/2018", ## dd/mm/yyyy

##  <#GAPDoc Label="PKGVERSIONDATA">
##  <!ENTITY VERSION "1.0.0">
##  <!ENTITY RELEASEDATE "05 December 2018">
##  <!ENTITY RELEASEYEAR "2018">
##  <#/GAPDoc>

PackageWWWHome :=
  Concatenation( "https://le-on.github.io/", ~.PackageName ),

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/le-on/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
SupportEmail := "junk@math.uni-sb.de",

ArchiveURL := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),

ArchiveFormats := ".tar.gz",


Persons := [
  rec(
    LastName      := "Junk",
    FirstNames    := "Luca Leon",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "junk@math.uni-sb.de",
    WWWHome       := "http://www.math.uni-sb.de/ag/weitze/",
    PostalAddress := Concatenation( [
                       "AG Weitze-Schmithüsen\n",
                       "FR 6.1 Mathematik\n",
                       "Universität des Saarlandes\n",
                       "D-66041 Saarbrücken" ] ),
    Place         := "Saarbrücken",
    Institution   := "Universität des Saarlandes"
  )

],

Status := "dev",

README_URL :=
  Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL :=
  Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),


AbstractHTML :=
  "This package provides a collection of algorithms for computing with \
  finite-index subgroups of (P)SL(2,Z).",


PackageDoc := rec(
  BookName  := "ModularGroup",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := ~.Subtitle,
),


Dependencies := rec(
  GAP := "4.5.3",

  NeededOtherPackages := [["GAPDoc", ">= 1.5"], ["CTblLib", ">= 1.2.2"]],

  SuggestedOtherPackages := [["Congruence", ">=1.1.1"]],

  ExternalConditions := []

),

AvailabilityTest := ReturnTrue,

BannerString := Concatenation(
    "----------------------------------------------------------------\n",
    "Loading  ModularGroup ", ~.Version, "\n",
    "by ",
    JoinStringsWithSeparator( List( Filtered( ~.Persons, r -> r.IsAuthor ),
                                    r -> Concatenation(
        r.FirstNames, " ", r.LastName, " (", r.WWWHome, ")\n" ) ), "   " ),
    "For help, type: ?ModularGroup package \n",
    "----------------------------------------------------------------\n" ),

TestFile := "tst/testall.g",

Keywords := ["PSL(2,Z)", "PSL2Z", "SL(2,Z)", "SL2Z", "modular group", "congruence subgroup"]

));
