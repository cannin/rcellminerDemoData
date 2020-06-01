context("Validation")

# TODO: Add better test
test_that("NCI-60 cell line names are consistent across molecular and drug data", {
  molDB <- rcellminer::getAllFeatureData(rcellminerData::molData)
  drugAct <- exprs(rcellminer::getAct(rcellminerData::drugData))
  
  expect_equal(TRUE, TRUE)
})
