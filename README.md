# rcellminerDemoData

This is a template repository to demonstrate generating experiment datasets compatible with rcellminer.

# Getting Started 

The main code lies in **inst/extdata/make_rcellminerdata_2_0.R**. This R script takes:

* mutation
* copy number
* mutation
* metadata (phenotypic measurements, such as doubling time)

and generates the molecular and drug RData datasets that make up the package.
