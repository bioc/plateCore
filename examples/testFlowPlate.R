######################################################################################################################
##
## Filename:testFlowPlate.R
##
## Author: straine
##
## Date: Oct 9, 2007
##
## Description:
##
######################################################################################################################

library(flowCore)
library(flowViz)
library(plateCore)

## Load the data from plateCore
data(plateCore)

ls()
#[1] "compensationSet"  "pbmcPlate"        "plateDescription" "wellAnnotation" 

## The pbmcPlate is a flowSet represent a large flow experiment ran on a 96 well
## plate.  The compensationSet is a flowSet containing the compensation controls.
## The plateDescription and wellAnnotation dataframes contain descriptions of the
## plate layout, and the contents of each well.

pbmcPlate
#A flowSet with 96 experiments.
#
#An object of class "AnnotatedDataFrame"
#rowNames: 0920609206.A01, 0920609206.A02, ..., 0920609206.H12  (96 total)
#varLabels and varMetadata description:
#		name: Name
#
#column names:
#		FSC-H SSC-H FL1-H FL2-H FL3-H FL1-A FL4-H Time


## The plateDescription name column must contain some mapping to the rowNames
## in the flowSet
plateDescription[1,]
#name Well.Id Row Column Well.Number
#1  A01     A01   A      1           1

## wellAnnotation describes what is in well
wellAnnotation[1,]
#Well.Id Channel       Dye Isotype.Group Sample.Type Target
#1     A01   FL1.H Alexa 488            12     Isotype     NA

## Create a flowPlate object 
platePBMC <- flowPlate(pbmcPlate,config=plateDescription,wellAnnot=wellAnnotation)
data <- flowPlate(pbmcPlate,config=plateDescription,wellAnnot=wellAnnotation)

## Calculate the compensation matrix using the spillover function from flowCore
comp.mat <- spillover(x=compensationSet,unstained=sampleNames(compensationSet)[5],patt=".*H",fsc="FSC.H",ssc="SSC.H",method="median")

## Since this data has already been compensated, set the off diagonal to 0
comp.mat[lower.tri(comp.mat)] <- 0
comp.mat[upper.tri(comp.mat)] <- 0

## Compensate the data.  Wells will only be compensated in channels for which they contain a dye.
platePBMC <- compensate(platePBMC,comp.mat)

## Truncate the values that are outside the detector range.

## Correct for autofluorescence

## QA/QC stuff

## Calculate the isotype gates.
platePBMC <- setContolGates(platePBMC,gateType="Isotype")

## Create a FlowJo workspace with the gates.



#
#wellAnnotation <- read.csv("wellAnnotationMP1.csv", as.is=TRUE,header=TRUE,stringsAsFactors=FALSE)
#x=flowPlate(pbmcPlate,config=plateDescription,wellAnnot=wellAnnotation)
#
#save(compensationSet,pbmcPlate,plateDescription,wellAnnotation,file="plateCore.rda")