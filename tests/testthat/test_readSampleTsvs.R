test_that("readSampleTsvs() works properly", {
    seqdataDir <- system.file("extdata/readSampleTsvs", package = "swissknife")
    
    ## Trying to read a non-existent file should give a warning and an empty data.frame
    expect_warning(df <- readSampleTsvs(seqdataDir = seqdataDir, 
                                        sampleIds = "readSampleTsvsEx4"))
    expect_equal(nrow(df), 0)
    expect_equal(ncol(df), 0)
    
    ## Check that the output is as expected for an existing file
    expect_equal(readSampleTsvs(seqdataDir = seqdataDir,
                                sampleIds = "readSampleTsvsEx1"),
                 data.frame(SAMPLE_ID = "readSampleTsvsEx1",
                            FIELD1 = "value1", FIELD2 = "value2",
                            TSV_FILE = "readSampleTsvsEx1.tsv",
                            stringsAsFactors = FALSE))
    
    ## Check that samples with multiple matches will not be read if keepMulti=FALSE
    expect_warning(df <- readSampleTsvs(seqdataDir = seqdataDir, 
                                        sampleIds = "readSampleTsvsEx",
                                        keepMulti = FALSE))
    expect_equal(nrow(df), 0)
    expect_equal(ncol(df), 0)
    
    ## Check that sample with multiple matches will be read properly if keepMulti=TRUE
    expect_warning(df <- readSampleTsvs(seqdataDir = seqdataDir, 
                                        sampleIds = "readSampleTsvsEx",
                                        keepMulti = TRUE))
    expect_equal(nrow(df), 3)
    expect_equal(ncol(df), 4)
    expect_equal(df$SAMPLE_ID, rep("readSampleTsvsEx", 3))
})