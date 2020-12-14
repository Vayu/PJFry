{ // rootlogon.C ...
  //
  // For ROOT, see: https://root.cern.ch/
  // The "rootlogon.C" macro, from the current working directory, is executed
  // automatically when the "root" starts (unless the "-n" option is used).
  //
  // std::cout << "... PJFry rootlogon.C ..." << std::endl;
  { // PJFry ...
    // try to find "libpjfry" (note: "LD_LIBRARY_PATH" may need to be set)
    char *p = gSystem->DynamicPathName("libpjfry", kTRUE);
    if (p && p[0]) { // "libpjfry" found ...
      if (gSystem->Load(p) >= 0) { // load the found "libpjfry" ...
        TString s(p);
        s.Remove(s.Last('/')); // strip "/libpjfry.*"
        // add PJFry include files' subdirectory
        if (s.EndsWith("/lib")) {
          // PJFry installed in "PREFIX/lib" and  "PREFIX/include"
          gInterpreter->AddIncludePath((s(0, s.Last('/')) + "/include").Data());
        } else if (s.EndsWith("/.libs")) {
          // PJFry in its source code directory tree
          gInterpreter->AddIncludePath((s(0, s.Last('/'))).Data());
        } else {
          // unknown directory structure (assume "in-place")
          gInterpreter->AddIncludePath(s.Data());
        }
        // any "local copies" of include files should take precedence
        // gInterpreter->AddIncludePath(".");
      } else { // ... "libpjfry" exists but it cannot be loaded ...
        std::cout << "... PJFry library cannot be loaded ..." << std::endl;
      } // ... load the found "libpjfry"
    } else { // ... "libpjfry" does not exist ...
      std::cout << "... PJFry library not found ..." << std::endl;
    } // ... "libpjfry" found
    delete p; // cleanup
  } // ... PJFry
} // ... rootlogon.C by Jacek M. Holeczek

