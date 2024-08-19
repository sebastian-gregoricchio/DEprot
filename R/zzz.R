.onAttach =
  function(libname, pkgname){
    msg = paste0("Loading required package: ", pkgname,"\n\n",
                 "If you are using DEprot in your work please cite:\n",
                 "'XYZ. XXJournal. 20XX.'")
    packageStartupMessage(msg)
  }
