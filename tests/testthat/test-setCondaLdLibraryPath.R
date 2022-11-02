test_that("setCondaLdLibraryPath works", {
    expect_error(setCondaLdLibraryPath(envPath = 1, condaLibPath = NULL),
                 "envPath must be a character scalar")
    expect_error(setCondaLdLibraryPath(envPath = c("path1", "path2"), 
                                       condaLibPath = NULL),
                 "envPath must be a character scalar")
    expect_error(setCondaLdLibraryPath(envPath = "path1", condaLibPath = 1),
                 "condaLibPath must be a character scalar or NULL")
    expect_error(setCondaLdLibraryPath(envPath = "path1", condaLibPath = "path"),
                 "The provided conda path doesn't exist")
    
    tDir <- tempdir()
    cPath <- tDir
    setCondaLdLibraryPath(envPath = tDir, condaLibPath = cPath)
    expect_true(file.exists(file.path(tDir, "etc", "conda", "activate.d", 
                                      "env_vars.sh")))
    rf <- readLines(file.path(tDir, "etc", "conda", "activate.d", 
                              "env_vars.sh"))
    expect_equal(length(rf), 4L)
    expect_equal(rf[1], "#!/bin/sh")
    expect_equal(rf[2], "")
    expect_equal(rf[3], "export TEMP_OLD_LD_LIBRARY_PATH=$LD_LIBRARY_PATH")
    expect_equal(rf[4], paste0("export LD_LIBRARY_PATH=", tDir, ":$LD_LIBRARY_PATH"))
    
    expect_true(file.exists(file.path(tDir, "etc", "conda", "deactivate.d", 
                                      "env_vars.sh")))
    rf <- readLines(file.path(tDir, "etc", "conda", "deactivate.d", 
                              "env_vars.sh"))
    expect_equal(length(rf), 4L)
    expect_equal(rf[1], "#!/bin/sh")
    expect_equal(rf[2], "")
    expect_equal(rf[3], "export LD_LIBRARY_PATH=$TEMP_OLD_LD_LIBRARY_PATH")
    expect_equal(rf[4], "unset TEMP_OLD_LD_LIBRARY_PATH")
})
