# overlayFilterResult-accessors.R
# 
# TODO: Add comment
#
# Author: straine
###############################################################################


setMethod("overlayFilterResult",signature("vector"),function(subSet,...) {
			temp <- new("overlayFilterResult")
			temp@subSet <- subSet
			return(temp)
		})


setMethod("overlay",signature("flowPlate"),function(data,type,channel,...) {
			if(type=="Negative.Control") {
				wellIds <- subset(data@wellAnnotation,Channel==channel)
				wellIds$eventCount <- sapply(wellIds$name,function(x) nrow(exprs(data@plateSet[[x]])))

				temp <- as(data@plateSet@frames,"list")				
				names(temp) <- sampleNames(data)
				
				frames <- lapply(names(temp),function(fileName) {
							if(fileName %in% wellIds$name) {
								negName <- wellIds[wellIds$name==fileName,c("Negative.Control","plateName")]
								if(length(negName)) negName <- wellIds[wellIds$Well.Id==negName[1,1] & wellIds$plateName==negName[1,2],"name"]
								if(length(negName) & negName %in% wellIds$name) {	
									exprs(temp[[fileName]]) <- rbind(exprs(temp[[fileName]]),exprs(temp[[negName]]))
								}
							}	
							temp[[fileName]]								
						})
				names(frames) <- names(temp)	
				
				plateSet <- as(frames,"flowSet")
				phenoData(plateSet) <- phenoData(data@plateSet)

				data@plateSet <- plateSet
				
				temp <- lapply(sampleNames(data@plateSet),function(sampName) {
							if(sampName %in% wellIds$name) {
								negName <- wellIds[wellIds$name==sampName,c("Negative.Control","plateName")]
								if(length(negName)) negName <- wellIds[wellIds$Well.Id==negName[1,1] & wellIds$plateName==negName[1,2],"name"]
								if(length(negName) & negName %in% wellIds$name) {
									return(overlayFilterResult(c(rep(0,wellIds[wellIds$name==sampName,"eventCount"]),
															rep(1,wellIds[wellIds$name==negName,"eventCount"]))))
								}
							}	
							overlayFilterResult(c(rep(0,nrow(exprs(data@plateSet[[sampName]])))))
						})
				
				names(temp) <- sampleNames(data@plateSet)
				
				data@overlayFilterResults <- temp
			}	
			data		
		})			

