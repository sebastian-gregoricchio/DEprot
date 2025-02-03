.onAttach =
  function(libname, pkgname){
    msg = paste0("Loading required package: ", pkgname,"\n\n",
                 "If you are using DEprot in your work please cite:\n",
                 "'XYZ. XXJournal. 20XX.'")
    packageStartupMessage(msg)


    # Check internet connection
    if (isFALSE(curl::has_internet())) {return(invisible(NULL))}

    # Store DEprot package versions
    check = suppressMessages(rvcheck::check_github("sebastian-gregoricchio/DEprot"))
    # When up-to-date == NA it means that it is a developmental version higher then the last release
    if (is.na(check$up_to_date)) {return(invisible(NULL))}

    # Check if DEprot is up-to-date and print a message if not up-to-date
    if (check$up_to_date == F) {
      message = paste("| The 'DEprot' package is not up-to-date. Installed version v",
                      check$installed_version, " --> v", check$latest_version,
                      " available.", sep="")
      command = '| To update DEprot type: devtools::install_github("sebastian-gregoricchio/DEprot", build_manual = T, build_vignettes = T)'

      message(paste("\033[0;37;44m",
                    paste(rep("-", nchar(command)+2), collapse = ""),
                    "\n", paste0(message, paste(rep(" ", nchar(command) - nchar(message)),  collapse = "")), " |\n",
                    command, " |",
                    "\n", paste(rep("-", max(c(nchar(message),nchar(command)))+2), collapse = ""),
                    "\033[0m",
                    sep = ""))
    }
  }
