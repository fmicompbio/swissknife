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
    
    ## Check that the output is as expected for multiple existing files
    expect_equal(readSampleTsvs(seqdataDir = seqdataDir,
                                sampleIds = c("readSampleTsvsEx1",
                                              "readSampleTsvsEx2")),
                 data.frame(SAMPLE_ID = c("readSampleTsvsEx1",
                                          "readSampleTsvsEx2"),
                            FIELD1 = c("value1", "value3"),
                            FIELD2 = c("value2", "value4"),
                            TSV_FILE = c("readSampleTsvsEx1.tsv",
                                         "readSampleTsvsEx2.tsv"),
                            stringsAsFactors = FALSE))
    
    ## Check that samples with multiple matches will not be kept if keepMulti=FALSE
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
    expect_equal(sort(df$FIELD1), sort(c("value1", "value3", "value5")))
})
