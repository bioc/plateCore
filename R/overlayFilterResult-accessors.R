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
				
				data@plateSet <- fsApply(data@plateSet,function(x) {
							sampName <- identifier(x)	
							if(sampName %in% wellIds$name) {
								negName <- wellIds[wellIds$name==sampName,"Negative.Control"]
								if(length(negName) && (negName %in% wellIds$Well.Id)) {	
									negName <- wellIds[wellIds$Well.Id==negName,"name"]
									exprs(x) <- rbind(exprs(x),exprs(data@plateSet[[negName]]))
								}
							}
							x	
						})
				
				temp <- lapply(sampleNames(data@plateSet),function(sampName) {
							if(sampName %in% wellIds$name) {
								negName <- wellIds[wellIds$name==sampName,"Negative.Control"]					
								if(length(negName) && (negName %in% wellIds$Well.Id)) {
									return(overlayFilterResult(c(rep(0,wellIds[wellIds$name==sampName,"eventCount"]),
															rep(1,wellIds[wellIds$Well.Id==negName,"eventCount"]))))
								}
							}	
							overlayFilterResult(c(rep(0,nrow(exprs(data@plateSet[[sampName]])))))
						})
				
				names(temp) <- sampleNames(data@plateSet)
				
				data@overlayFilterResults <- temp
			}	
			data		
		})			

