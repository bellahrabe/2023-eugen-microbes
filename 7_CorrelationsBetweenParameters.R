basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
package.list <- setdiff(package.list,basic.packages)
if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
rm(list=ls()); closeAllConnections(); gc()

library(dplyr)

# # to create plots colored and grouped as in VisuaR analysis ###
VisuaRProjectName <- ""
PathToVisuaRAnalysis <- file.path("")

filename <- "water-t-vs-dna-6max" # name of analysis. e.g., temperature-vs-depth

source(file.path(PathToVisuaRAnalysis,VisuaRProjectName,paste0(VisuaRProjectName,"_UserInput.R",sep="")))

ylimitsset = "N"
ylims.plot <- c(NULL,NULL) # set to c(NULL,NULL) if no y axis rescaling is wanted or e.g. c(0,25) 

StartParamX <- "N" # set to 0 if you want to start the scale at 0 (1 for log scale), or to "N" if automatically start at the lowest value
StartParamLogX <- "N" # set to same as StartParamX or to 1 if StartParamX is set to 0
EndParamX <- "N" # set end of x axis scale

colordots <- "VisuaR" # set to VisuaR if dots in scatterplot should be colored only by the VisuaR grouping, give continuous variable to color by that
my_colors <- colorspace::sequential_hcl(n=7,palette="Plasma") # colorscale used for the continuous variable

# # additional metadata (than what you used in your VisuaR analysis)? should be in csv format and sample names as first column, parameters in other columns
# # after the column name the first three rows in the contextual data should be 
# # 1. parameter name: a name without special characters, this name will be used in the plot naming
# # 2. legend name: this name will be plotted as legend header 
# # 3. parameter type: continuous or discrete 
maincontext <- data.table::fread(input="path/to/data.csv")

# # read in from VisuaR output
diversityindices <- data.table::fread(input=paste0(PathToVisuaRAnalysis,"/",VisuaRProjectName,"/Alpha_Diversity/",VisuaRProjectName,"_DiversityIndices.csv"))
# # add first three rows for plotting
col_types <- sapply(diversityindices, function(col) {
  if (all(is.numeric(col))) {
    "continuous"
  } else {
    "discrete"
  }
})
col_types[1] <- "Dummy3"
add_rows <- data.frame(Dummy1 = colnames(diversityindices), Dummy2 = colnames(diversityindices),
                       Dummy3= col_types)
add_rows[1,1] <- "Dummy1"
add_rows[1,2] <- "Dummy2"
diversityindices <- rbind(
  t(add_rows),
  diversityindices
)
rm(add_rows)

# Replace spaces, brackets, hyphens, and % in the first row with nothing
diversityindices[1, ] <- lapply(diversityindices[1, ], function(x) gsub("[\\(\\)_% -]", "", x))

# # reads in used contextual data for ViusaR analysis togerther with the calculated alpha diversity values
# # diversity indices (above) contains some additional info though and nicer naming
t.M.diversity.ord <- readRDS(file=file.path(PathToVisuaRAnalysis,VisuaRProjectName,"Alpha_Diversity",paste(VisuaRProjectName,"_DiversityIndices-t.rds",sep="")))

# # just in case some variables where added
t.M.diversity.ord <- dplyr::left_join(t.M.diversity.ord,diversityindices,by=c("SampleNameVisuaR"="V1"))

t.M.diversity.ord <- dplyr::left_join(t.M.diversity.ord,maincontext,by=c("SampleNameVisuaR"="SampleNameVisuaR"))
rm(diversityindices,maincontext)

ParameterToPlotY <- "" # e.g. "SST"
ParameterToPlotX <- "" # e.g. "DNA-DNA-concentration-in-seawater-ng-per-ml"

excludeSamples <-  "Y" # if set to yes you can define a category you want to exclude, e.g. "soil" from your plots. 
CategoryToExclude <- "air" #Optional. this is a categories in your metadata file. If you do not want to exclude a category leave it empty
FoundInColumnName <- "VisuaR-SampleType" #Optional. the metadata file's column name where the 'CategoryToExclude' variable can be found. Leave empty when not needed.

# # do not change anything from here on ###################
# # data read-in and modification ####
# # read in contextual data (last column is color vector as used in the VisuaR analysis)

# # exclude particular samples
if (excludeSamples == "Y") {
  t.M.diversity.ord=t.M.diversity.ord[grep(CategoryToExclude,t.M.diversity.ord[,which(names(t.M.diversity.ord)==FoundInColumnName)],invert=T),] # excludes samples belonging to the respective category
}

# # create folder for plot output
PathToParameterPlots=file.path(PathToVisuaRAnalysis,'parameterplots') # creates directories for the current VisuaR analysis
dir.create(PathToParameterPlots)
PathToParameterPlots=file.path(PathToVisuaRAnalysis,'parameterplots',filename) # creates directories for the current VisuaR analysis
dir.create(PathToParameterPlots)

# # prepare datatable for plotting
# # data table needs to be with names saved as factor in order to keep the order of samples while plotting
t.M.diversity.ord.fac<-as.data.frame(t.M.diversity.ord)
t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ColumnForAnalysisNames)] <- factor(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ColumnForAnalysisNames)], levels=t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ColumnForAnalysisNames)])

# # subset so noNAs
t.M.diversity.ord.fac <- t.M.diversity.ord.fac[!is.na(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ParameterToPlotY)]),]
t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ParameterToPlotY)] <- as.numeric(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ParameterToPlotY)])

# # subset so noNAs
t.M.diversity.ord.fac <- t.M.diversity.ord.fac[!is.na(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ParameterToPlotX)]),]
t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ParameterToPlotX)] <- as.numeric(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==ParameterToPlotX)])
if (colordots != "VisuaR") {
  t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==colordots)] <- as.numeric(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==colordots)])
}

yparamrange <- (max(t.M.diversity.ord.fac[, which(names(t.M.diversity.ord.fac) == ParameterToPlotY)])-min(t.M.diversity.ord.fac[, which(names(t.M.diversity.ord.fac) == ParameterToPlotY)]))

# calculates the number of pairwise comparisons for bonferroni correction (0 when 1 group, 1 when 2 groups,3 when 3 groups, 6 when 4 groups etc.)
numpairs<-length(unique(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==Grouping1)]))*(length(unique(t.M.diversity.ord.fac[,which(names(t.M.diversity.ord.fac)==Grouping1)]))-1)/2

# # Parameter vs Parameter Plots #####
t.M.diversity.ord.fac.sub=subset(t.M.diversity.ord.fac,!is.na(t.M.diversity.ord.fac[,ParameterToPlotY])&!is.na(t.M.diversity.ord.fac[,ParameterToPlotX])) # creates a Sample by Metadata file of samples which have NAs in one of the selected Groupings
p1 <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, 
                      ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)], 
                                   y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)]))+
  ggplot2::geom_point(fill=as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]), size=2.5, stroke=0.5, shape=21, color="black") +
  ggplot2::labs(x = ParameterToPlotX, y = ParameterToPlotY)+
  ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
  ggplot2::geom_smooth(method='lm',formula = y~x,linetype='solid',color='black',fill="grey",se=T,fullrange=T,level=0.95)+ #smoothing methods are auto, lm, glm, gam, loess, rlm. loess: default value for small number of observations. lm fits are linear model, se=T:displays confidence interval around smooth
  ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
  ggpubr::stat_cor(label.x.npc = 0.7,label.y.npc = 0.9,size=3) + # note: his shows Pearsons Rho and not R2 #label.x = 4,label.y = 13
  ggpubr::stat_regline_equation(label.x.npc = 0.7,label.y.npc = 0.85,size=3)+ #label.x=4,label.y=11
  ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                 axis.text.y=ggplot2::element_text(colour="black",size=10),
                 axis.text.x=ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                 panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                 panel.grid.minor.y = ggplot2::element_line(),
                 panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                 axis.title = ggplot2::element_text(colour="black",size=12),
                 axis.line= ggplot2::element_line(colour="black"),
                 panel.background = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                 legend.title = ggplot2::element_blank(),
                 legend.box.background = ggplot2::element_rect(colour="black"),
                 legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                 legend.key= ggplot2::element_rect(fill="white"),
                 legend.text= ggplot2::element_text(size=14),
                 legend.position = 'none',
                 axis.ticks= ggplot2::element_line(colour = "black"),
                 plot.background = ggplot2::element_rect(size=c(1,1)))
if (StartParamX =="N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = ylims.plot)
  } else {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])))
  }
} else if (StartParamX != "N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = ylims.plot,
                                        xlim = c(StartParamX,
                                                 NA))
  } else {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                        xlim = c(StartParamX,
                                                 NA))
  }
} else if (StartParamX != "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = ylims.plot,
                                        xlim = c(StartParamX,
                                                 EndParamX))
  } else {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                        xlim = c(StartParamX,
                                                 EndParamX))
  }
} else if (StartParamX == "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = ylims.plot,
                                        xlim = c(NA,
                                                 EndParamX))
  } else {
    p1 <- p1 + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                        xlim = c(NA,
                                                 EndParamX))
  }
}
p1
  
if (colordots != "VisuaR") {
  p1filled <- ggplot2::ggplot(t.M.diversity.ord.fac.sub,
                              ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)], 
                                           y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)],
                                           fill=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==colordots)]))+
    ggplot2::geom_point(size=2.5, stroke=0.5, shape=21, color="black") +
    ggplot2::scale_fill_gradientn(colors=my_colors) +
    ggplot2::labs(x = ParameterToPlotX, y = ParameterToPlotY, fill= colordots)+
    ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                     axis.text.y=ggplot2::element_text(colour="black",size=10),
                     axis.text.x=ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                     panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                     panel.grid.minor.y = ggplot2::element_line(),
                     panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                     axis.title = ggplot2::element_text(colour="black",size=12),
                     axis.line= ggplot2::element_line(colour="black"),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                     legend.title = ggplot2::element_text(size=8),
                     legend.box.background = ggplot2::element_blank(),
                     legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                     legend.key= ggplot2::element_rect(fill="white"),
                     legend.text= ggplot2::element_text(size=8),
                     axis.ticks= ggplot2::element_line(colour = "black"),
                     plot.background = ggplot2::element_rect(size=c(1,1)))
    
  if (StartParamX =="N" & EndParamX == "N") {
    if (ylimitsset == "Y") {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = ylims.plot)
    } else {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                               max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])))
    }
  } else if (StartParamX != "N" & EndParamX == "N") {
    if (ylimitsset == "Y") {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = ylims.plot,
                                                      xlim = c(StartParamX,
                                                               NA))
    } else {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                               max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                                      xlim = c(StartParamX,
                                                               NA))
    }
  } else if (StartParamX != "N" & EndParamX != "N") {
    if (ylimitsset == "Y") {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = ylims.plot,
                                                      xlim = c(StartParamX,
                                                               EndParamX))
    } else {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                               max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                                      xlim = c(StartParamX,
                                                               EndParamX))
    }
  } else if (StartParamX == "N" & EndParamX != "N") {
    if (ylimitsset == "Y") {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = ylims.plot,
                                                      xlim = c(NA,
                                                               EndParamX))
    } else {
      p1filled <- p1filled + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                               max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                                      xlim = c(NA,
                                                               EndParamX))
    }
  }
  pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-by-",colordots,".pdf",sep=""),sep = '')),height=4,width=5.2,useDingbats=F); print(p1filled); dev.off()
}
  
p1pseudolog <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)], y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)]))+
    ggplot2::geom_point(fill=as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]), size=2.5, stroke=0.5, shape=21, color="black")+
    ggplot2::labs(x = ParameterToPlotX, y = ParameterToPlotY)+
    ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_x_continuous(trans='pseudo_log',breaks = c(0,1,2,3,4,10,100,500,1000,10000,100000,500000,1000000,5000000,10000000,100000000,1000000000),labels = scales::label_number_auto()) +
    ggpubr::stat_cor(label.x.npc = 0.7,label.y.npc = 0.9,size=3) + # note: his shows Pearsons Rho and not R2 #label.x = 4,label.y = 13
    ggpubr::stat_regline_equation(label.x.npc = 0.7,label.y.npc = 0.85,size=3)+ #label.x=4,label.y=11
    ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                   axis.text.y= ggplot2::element_text(colour="black",size=10),
                   axis.text.x= ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                   panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                   panel.grid.minor.y = ggplot2::element_line(),
                   panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                   axis.title = ggplot2::element_text(colour="black",size=12),
                   axis.line= ggplot2::element_line(colour="black"),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                   legend.title = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_rect(colour="black"),
                   legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                   legend.key= ggplot2::element_rect(fill="white"),
                   legend.text= ggplot2::element_text(size=14),
                   legend.position = 'none',
                   axis.ticks= ggplot2::element_line(colour = "black"),
                   plot.background = ggplot2::element_rect(size=c(1,1)))
  
if (StartParamX =="N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = ylims.plot)
  } else {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])))
  }
} else if (StartParamX != "N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = ylims.plot,
                                                          xlim = c(StartParamX,
                                                                   NA))
  } else {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                                          xlim = c(StartParamX,
                                                                   NA))
  }
} else if (StartParamX != "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = ylims.plot,
                                                          xlim = c(StartParamX,
                                                                   EndParamX))
  } else {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                                          xlim = c(StartParamX,
                                                                   EndParamX))
  }
} else if (StartParamX == "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = ylims.plot,
                                                          xlim = c(NA,
                                                                   EndParamX))
  } else {
    p1pseudolog <- p1pseudolog + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                                          xlim = c(NA,
                                                                   EndParamX))
  }
}
p1pseudolog
  
p1log <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)], y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)]))+
    ggplot2::geom_point(fill=as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]), size=2.5, stroke=0.5, shape=21, color="black")+
    ggplot2::labs(x = ParameterToPlotX, y = ParameterToPlotY)+
    ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_x_continuous(trans='log10',breaks = c(1,2,3,4,10,100,500,1000,10000,100000,500000,1000000,5000000,10000000,100000000,1000000000),labels = scales::label_number_auto()) +
    ggpubr::stat_cor(label.x.npc = 0.7,label.y.npc = 0.9,size=3) + # note: his shows Pearsons Rho and not R2 #label.x = 4,label.y = 13
    ggpubr::stat_regline_equation(label.x.npc = 0.7,label.y.npc = 0.85,size=3)+ #label.x=4,label.y=11
    ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                   axis.text.y= ggplot2::element_text(colour="black",size=10),
                   axis.text.x= ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                   panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                   panel.grid.minor.y = ggplot2::element_line(),
                   panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                   axis.title = ggplot2::element_text(colour="black",size=12),
                   axis.line= ggplot2::element_line(colour="black"),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                   legend.title = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_rect(colour="black"),
                   legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                   legend.key= ggplot2::element_rect(fill="white"),
                   legend.text= ggplot2::element_text(size=14),
                   legend.position = 'none',
                   axis.ticks= ggplot2::element_line(colour = "black"),
                   plot.background = ggplot2::element_rect(size=c(1,1)))
if (StartParamX =="N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = ylims.plot)
  } else {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                       max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])))
  }
} else if (StartParamX != "N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = ylims.plot,
                                              xlim = c(StartParamX,
                                                       NA))
  } else {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                       max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                              xlim = c(StartParamX,
                                                       NA))
  }
} else if (StartParamX != "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = ylims.plot,
                                              xlim = c(StartParamX,
                                                       EndParamX))
  } else {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                       max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                              xlim = c(StartParamX,
                                                       EndParamX))
  }
} else if (StartParamX == "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = ylims.plot,
                                              xlim = c(NA,
                                                       EndParamX))
  } else {
    p1log <- p1log + ggplot2::coord_cartesian(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)]),
                                                       max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])),
                                              xlim = c(NA,
                                                       EndParamX))
  }
}
p1log
  
# # create plots with bars coming from the left instead of dots
p1bar <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)], x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)]))+
    ggplot2::geom_bar(stat="identity",fill=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")],width=yparamrange*0.0285714,alpha=0.8,color="black")+
    ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                   axis.text.y=ggplot2::element_text(colour="black",size=10),
                   axis.text.x=ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                   panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                   panel.grid.minor.y = ggplot2::element_line(),
                   panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                   axis.title = ggplot2::element_text(colour="black",size=12),
                   axis.line= ggplot2::element_line(colour="black"),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                   legend.title = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_rect(colour="black"),
                   legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                   legend.key= ggplot2::element_rect(fill="white"),
                   legend.text= ggplot2::element_text(size=14),
                   legend.position = 'none',
                   axis.ticks= ggplot2::element_line(colour = "black"),
                   plot.background = ggplot2::element_rect(size=c(1,1)))
  
if (StartParamX =="N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1bar <- p1bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1bar <- p1bar + 
      ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX != "N" & EndParamX == "N") {
  if (ylimitsset =="Y") {
    p1bar <- p1bar + 
      ggplot2::coord_flip(xlim =ylims.plot,
                          ylim = c(StartParamX,
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1bar <- p1bar + 
      ggplot2::coord_flip(ylim = c(StartParamX,
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX != "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1bar <- p1bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(StartParamX,
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1bar <- p1bar + 
      ggplot2::coord_flip(ylim = c(StartParamX,
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX == "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1bar <- p1bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1bar <- p1bar + 
      ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
}
p1bar

p1pseudolog_bar <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)], y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)]))+
  ggplot2::geom_bar(stat="identity",fill=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")],width=yparamrange*0.0285714,,alpha=0.8,color="black")+
  ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
  ggplot2::scale_y_continuous(trans='pseudo_log',breaks = c(0,1,2,3,4,10,100,500,1000,10000,100000,500000,1000000,5000000,10000000,100000000,1000000000),labels = scales::label_number_auto()) +
  ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                 axis.text.y= ggplot2::element_text(colour="black",size=10),
                 axis.text.x= ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                 panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                 panel.grid.minor.y = ggplot2::element_line(),
                 panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                 axis.title = ggplot2::element_text(colour="black",size=12),
                 axis.line= ggplot2::element_line(colour="black"),
                 panel.background = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                 legend.title = ggplot2::element_blank(),
                 legend.box.background = ggplot2::element_rect(colour="black"),
                 legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                 legend.key= ggplot2::element_rect(fill="white"),
                 legend.text= ggplot2::element_text(size=14),
                 legend.position = 'none',
                 axis.ticks= ggplot2::element_line(colour = "black"),
                 plot.background = ggplot2::element_rect(size=c(1,1)))

if (StartParamX =="N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX != "N" & EndParamX == "N") {
  if (ylimitsset =="Y") {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(xlim =ylims.plot,
                          ylim = c(StartParamX,
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(ylim = c(StartParamX,
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX != "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(StartParamX,
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(ylim = c(StartParamX,
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX == "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1pseudolog_bar <- p1pseudolog_bar + 
      ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
}
p1pseudolog_bar

p1log_bar <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)], y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)]))+
    ggplot2::geom_bar(stat="identity",fill=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")],width=yparamrange*0.0285714,,alpha=0.8,color="black")+
    ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_y_continuous(trans='log10',breaks = c(1,2,3,4,10,100,500,1000,10000,100000,500000,1000000,5000000,10000000,100000000,1000000000),labels = scales::label_number_auto()) +
    ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                   axis.text.y= ggplot2::element_text(colour="black",size=10),
                   axis.text.x= ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                   panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                   panel.grid.minor.y = ggplot2::element_line(),
                   panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
                   axis.title = ggplot2::element_text(colour="black",size=12),
                   axis.line= ggplot2::element_line(colour="black"),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                   legend.title = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_rect(colour="black"),
                   legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                   legend.key= ggplot2::element_rect(fill="white"),
                   legend.text= ggplot2::element_text(size=14),
                   legend.position = 'none',
                   axis.ticks= ggplot2::element_line(colour = "black"),
                   plot.background = ggplot2::element_rect(size=c(1,1)))
if (StartParamX =="N" & EndParamX == "N") {
  if (ylimitsset == "Y") {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX != "N" & EndParamX == "N") {
  if (ylimitsset =="Y") {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(xlim =ylims.plot,
                          ylim = c(StartParamX,
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(ylim = c(StartParamX,
                                   max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX != "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(StartParamX,
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(ylim = c(StartParamX,
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
} else if (StartParamX == "N" & EndParamX != "N") {
  if (ylimitsset == "Y") {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(xlim = ylims.plot,
                          ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  } else {
    p1log_bar <- p1log_bar + 
      ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                   EndParamX)) +
      ggplot2::labs(x = ParameterToPlotY, y = ParameterToPlotX)
  }
}

# # calculate Rho and p values for the correlation of ParamterToPlot2 and ParameterToPlotY (NA for groups with less than 3 samples)
cor.values <- as.data.frame(matrix(nrow = length(unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")])), ncol = 2))
colnames(cor.values) <- c("Rho","p")

cor.colors <- vector(length=length(unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")])))

for (i in 1:length(unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")]))) {
    t.M.diversity.ord.fac.sub1<-t.M.diversity.ord.fac.sub[t.M.diversity.ord.fac.sub[["M.colvec.ord"]] == unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")])[[i]],c( which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX), which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY),which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord"),which(names(t.M.diversity.ord.fac.sub) == Grouping1)) ]
    if (nrow(t.M.diversity.ord.fac.sub1)<3) {
      rownames(cor.values)[i] <-  t.M.diversity.ord.fac.sub1[1,4]
      cor.colors[i] <- t.M.diversity.ord.fac.sub1[1,3]
    } else if (nrow(t.M.diversity.ord.fac.sub1)>2) {
      rownames(cor.values)[i] <-  t.M.diversity.ord.fac.sub1[1,4]
      cor.values[i,1] <- as.numeric(cor.test(t.M.diversity.ord.fac.sub1[,1],t.M.diversity.ord.fac.sub1[,2]) [4]) # extracts Rho
      cor.values[i,2] <- as.numeric(cor.test(t.M.diversity.ord.fac.sub1[,1],t.M.diversity.ord.fac.sub1[,2]) [3]) # extracts p
      cor.colors[i] <- t.M.diversity.ord.fac.sub1[1,3]
    }
}
data.table::fwrite(cor.values,file=file.path(PathToParameterPlots,paste(paste0(filename,"_corvalues.csv",sep=""),sep = '')),col.names = T,row.names = T)
rm(t.M.diversity.ord.fac.sub1)
  
# # create color vector for regline: only keeps groups which occur at least twice
reg.colors <- t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")][duplicated(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")])|duplicated(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")],fromLast = T)]
reg.colors.unique <- unique (reg.colors)
# # get grouped regline equations and R-squared
reg.values <- as.data.frame(matrix(nrow = length(unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")])), ncol = 2))
colnames(reg.values) <- c("Intercept","x")

for (i in 1:length(unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")]))) {
  t.M.diversity.ord.fac.sub1<-t.M.diversity.ord.fac.sub[t.M.diversity.ord.fac.sub[["M.colvec.ord"]] == unique(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == "M.colvec.ord")])[[i]],c( which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX),which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY),which(names(t.M.diversity.ord.fac.sub) == Grouping1)) ]
  colnames(t.M.diversity.ord.fac.sub1) <- c("x","y","Grouping1")
  if (nrow(t.M.diversity.ord.fac.sub1)<2) {
    rownames(reg.values)[i] <- t.M.diversity.ord.fac.sub1[1,3]
  } else if (nrow(t.M.diversity.ord.fac.sub1)>1) {
    reg.values[i,1] <- lm(y~x,t.M.diversity.ord.fac.sub1)[[1]][1] # returns the y intercept
    reg.values[i,2] <- lm(y~x,t.M.diversity.ord.fac.sub1)[[1]][2] # returns the x value
    rownames(reg.values)[i] <- t.M.diversity.ord.fac.sub1[1,3]
  }
}
data.table::fwrite(reg.values,file=file.path(PathToParameterPlots,paste(paste0(filename,"_regvalues.csv",sep=""),sep = '')),col.names = T,row.names = T)
rm(t.M.diversity.ord.fac.sub1)
  
# # plots Parameter 2 (on x axis as boxplots so it can be plotted on top of scatterplot)
p2 <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)],x=t.M.diversity.ord.fac.sub[,Grouping1])) + 
    ggplot2::geom_boxplot(fill=rev(unique(as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]))))  
if (StartParamX =="N" & EndParamX == "N") {
  p2 <- p2 +
    ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) 
} else if (StartParamX != "N" & EndParamX == "N") {
  p2 <- p2 +
    ggplot2::coord_flip(ylim = c(StartParamX,
                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)])))
} else if (StartParamX != "N" & EndParamX != "N") {
  p2 <- p2 +
    ggplot2::coord_flip(ylim = c(StartParamX,
                                 EndParamX))
} else if (StartParamX == "N" & EndParamX != "N") {
  p2 <- p2 +
    ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                 EndParamX))
}
  
p2 <- p2 +
    ggplot2::scale_x_discrete(limits=c(rev(unique(t.M.diversity.ord.fac.sub[,Grouping1]))),labels=NULL)+
    ggplot2::xlab(NULL)+ # to have blank axes
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line= ggplot2::element_line(colour="black"),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
      axis.ticks.x= ggplot2::element_line(colour = "black"),
      axis.ticks.y= ggplot2::element_blank(),
      legend.position = "none",
      plot.background = ggplot2::element_rect(size=c(1,1)))
print(p2)
  
# # plots Parameter 2 (on x axis as boxplots so it can be plotted on top of scatterplot) with log10 scale
p2log <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)],x=t.M.diversity.ord.fac.sub[,Grouping1])) + 
  ggplot2::geom_boxplot(fill=rev(unique(as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")])))) 
if (StartParamX =="N" & EndParamX == "N") {
  p2log <- p2log +
    ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) 
} else if (StartParamX != "N" & EndParamX == "N") {
  p2log <- p2log + 
    ggplot2::coord_flip(ylim = c(StartParamLogX,
                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)])))
} else if (StartParamX != "N" & EndParamX != "N") {
  p2log <- p2log + 
    ggplot2::coord_flip(ylim = c(StartParamLogX,
                                 EndParamX))
} else if (StartParamX == "N" & EndParamX != "N") {
  p2log <- p2log +
    ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                 EndParamX)) 
}
p2log <- p2log +
  ggplot2::scale_x_discrete(limits=c(rev(unique(t.M.diversity.ord.fac.sub[,Grouping1]))),labels=NULL)+
  ggplot2::scale_y_continuous(trans='log10') +
  ggplot2::xlab(NULL)+ # to have blank axes
  ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line= ggplot2::element_line(colour="black"),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
      axis.ticks.x = ggplot2::element_line(colour = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.background = ggplot2::element_rect(size=c(1,1)))
# # plots Parameter 2 (on x axis as boxplots so it can be plotted on top of scatterplot) with pseudo-log10 scale
p2pseudolog <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)],x=t.M.diversity.ord.fac.sub[,Grouping1])) + 
  ggplot2::geom_boxplot(fill=rev(unique(as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]))))  
if (StartParamX =="N" & EndParamX == "N") {
  p2pseudolog <- p2pseudolog +
    ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]))) 
} else if (StartParamX != "N" & EndParamX == "N") {
  p2pseudolog <- p2pseudolog +
    ggplot2::coord_flip(ylim = c(StartParamX,
                                 max(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)])))
} else if (StartParamX != "N" & EndParamX != "N") {
  p2pseudolog <- p2pseudolog +
    ggplot2::coord_flip(ylim = c(StartParamX,
                                 EndParamX))
} else if (StartParamX == "N" & EndParamX != "N") {
  p2pseudolog <- p2pseudolog +
    ggplot2::coord_flip(ylim = c(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotX)]),
                                 EndParamX)) 
}
p2pseudolog <- p2pseudolog +
    ggplot2::scale_x_discrete(limits=c(rev(unique(t.M.diversity.ord.fac.sub[,Grouping1]))),labels=NULL)+
    ggplot2::scale_y_continuous(trans='pseudo_log',breaks = c(1,10,100,1000,10000,100000,1000000),labels = scales::label_number_auto()) +
    ggplot2::xlab(NULL)+ # to have blank axes
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour="grey",size=0.25),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line= ggplot2::element_line(colour="black"),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
      axis.ticks.x = ggplot2::element_line(colour = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none",
      plot.background = ggplot2::element_rect(size=c(1,1)))
print(p2pseudolog)
  
if(numpairs>=1) {
    p2n <- p2 + 
      ggsignif::geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==Grouping1)]))),2,simplify=F, test="wilcox.test"), # default: wilcox.test()
                            vjust = 0.6,
                            step_increase=0.06,
                            size=0.2,
                            textsize=3,
                            tip_length=0.01,
                            map_signif_level=c("***"=0.001/numpairs,"**"=0.01/numpairs,"*"=0.05/numpairs," "=2))+
      ggplot2::coord_flip()
    print(p2n) 
    p2nuncor <- p2 + 
      ggsignif::geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==Grouping1)]))),2,simplify=F, test="wilcox.test"), # default: wilcox.test()
                            vjust = 0.6,
                            step_increase=0.06,
                            size=0.2,
                            textsize=3,
                            tip_length=0.01,
                            map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2)) +
      ggplot2::coord_flip()
    print(p2nuncor) 
}

p3 <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)],x=t.M.diversity.ord.fac.sub[,Grouping1])) + 
    ggplot2::geom_boxplot(fill=unique(as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]))) + 
    ggplot2::xlab(ParameterToPlotY)+
    ggplot2::ylab("")+
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.text.x= ggplot2::element_text(angle=0,colour="white",size=10),#angle=45,hjust=1
                   axis.text.y= ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(colour="white",size=12),
                   axis.title.y = ggplot2::element_blank(),
                   axis.line= ggplot2::element_line(colour="black"),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                   panel.grid.minor.y = ggplot2::element_line(),
                   panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                   legend.box.background = ggplot2::element_rect(colour="black"),
                   legend.box.margin = ggplot2::margin(0.18,0.1,0.1,0.5),
                   legend.key = ggplot2::element_rect(fill="white"),
                   legend.text= ggplot2::element_text(size=14),
                   legend.position = 'none',
                   axis.ticks= ggplot2::element_line(colour = "black"),
                   axis.ticks.x= ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(size=c(1,1)))
 if (ylimitsset == "Y") {
   p3 <- p3 + ggplot2::coord_cartesian(ylim = ylims.plot)
 } 
 print(p3)
  
 if (numpairs>=1) {
    p3n <- p3 + EnvStats::stat_n_text(y.expand.factor = 0,y.pos =((ggplot2::layer_scales(p1)$y$range$range[1])-(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])))) + 
      ggsignif::geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==Grouping1)]))),2,simplify=F, test="wilcox.test"), # default: wilcox.test()
                            vjust = 0.6,
                            step_increase=0.06,
                            size=0.2,
                            textsize=3,
                            tip_length=0.01,
                            map_signif_level=c("***"=0.001/numpairs,"**"=0.01/numpairs,"*"=0.05/numpairs," "=2)) 
    
    print(p3n)  
    p3nuncor <- p3 + EnvStats::stat_n_text(y.expand.factor = 0,y.pos =((ggplot2::layer_scales(p1)$y$range$range[1])-(min(t.M.diversity.ord.fac.sub[, which(names(t.M.diversity.ord.fac.sub) == ParameterToPlotY)])))) + 
      ggsignif::geom_signif(comparisons=combn(sort(as.character(unique(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==Grouping1)]))),2,simplify=F, test="wilcox.test"), # default: wilcox.test()
                            vjust = 0.6,
                            step_increase=0.06,
                            size=0.2,
                            textsize=3,
                            tip_length=0.01,
                            map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05," "=2)) 
    print(p3nuncor)  
 }
  
p4 <- ggplot2::ggplot(t.M.diversity.ord.fac.sub, ggplot2::aes(x=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotX)], y=t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)==ParameterToPlotY)])) +
    ggplot2::geom_point(ggplot2::aes(color=as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")])),size=2.5, stroke=0.5, shape=15) +
    ggplot2::scale_color_identity(name=colnames(t.M.diversity.ord.fac.sub)[which(colnames(t.M.diversity.ord.fac.sub)==Grouping1)],
                                  breaks = c(unique(as.character(t.M.diversity.ord.fac.sub[,which(names(t.M.diversity.ord.fac.sub)=="M.colvec.ord")]))),
                                  labels = c(unique(t.M.diversity.ord.fac.sub[,which(colnames(t.M.diversity.ord.fac.sub)==Grouping1)])),
                                  guide = "legend")+
    ggplot2::geom_smooth(method='lm',formula = y~x,linetype='solid',color='black',fill="grey",se=T,fullrange=T,level=0.95)+ #smoothing methods are auto, lm, glm, gam, loess, rlm. loess: default value for small number of observations. lm fits are linear model, se=T:displays confidence interval around smooth
    ggplot2::scale_y_continuous(breaks = c(0,25,50,75,100))+
    ggplot2::scale_x_continuous(labels = scales::label_number_auto())+
    ggpubr::stat_cor(label.x = 40,label.y = 80) + # note: his shows Pearsons Rho and not R2
    ggpubr::stat_regline_equation(label.x=40,label.y=90)+
    ggplot2::coord_trans(ylim = c(0,100)) + 
    ggplot2::theme(axis.text = ggplot2::element_text(colour="black"),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   axis.text.x= ggplot2::element_text(angle=0,colour="black",size=10),#angle=45,hjust=1
                   axis.text.y= ggplot2::element_text(colour="black",size=10),
                   axis.title.y = ggplot2::element_text(colour="black",size=12),
                   axis.title.x = ggplot2::element_blank(),
                   axis.line= ggplot2::element_line(colour="black"),
                   panel.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour="grey",size=0.25),
                   panel.grid.minor.y = ggplot2::element_line(),
                   panel.border = ggplot2::element_rect(fill=NA,colour = "black",size=0.5),
                   legend.box.background = ggplot2::element_rect(colour="white"),
                   legend.box.margin = ggplot2::margin(0.1,0.1,0.1,0.1),
                   legend.key= ggplot2::element_rect(fill="white"),
                   legend.text= ggplot2::element_text(size=10),
                   legend.position = 'right',
                   axis.ticks= ggplot2::element_line(colour = "black"),
                   plot.background = ggplot2::element_rect(size=c(1,1)))
p4
legend <- cowplot::get_legend(p4)
  
# # Modify margin c(top, right, bottom, left) to reduce the distance between plots
# #and align G1 density with the scatterplot
p2 = p2 + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.1, 0, 2), "cm"))
if (numpairs>=1) {
    p2n = p2n + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.1, 0, 2), "cm"))
    p2nuncor <- p2nuncor + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.1, 0, 2), "cm"))
}
p2log <- p2log + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.1, 0, 2), "cm"))
p2pseudolog <- p2pseudolog + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 0.1, 0, 2), "cm"))

p1 = p1 + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.5, 0.5), "cm"))
p1log <- p1log + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.5, 0.5), "cm"))
p1pseudolog <- p1pseudolog + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.5, 0.5), "cm"))
p1bar <- p1bar + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.5, 0.5), "cm"))
p1log_bar <- p1log_bar + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.5, 0.5), "cm"))
p1pseudolog_bar <- p1pseudolog_bar + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.5, 0.5), "cm"))
  
p3 = p3 + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.5, 0.5, 0), "cm"))
if (numpairs>=1) {
  p3n = p3n + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.5, 0.5, 0), "cm"))
  p3nuncor <- p3nuncor + ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.5, 0.5, 0), "cm"))
}
  
# Combine all plots together and crush graph density with rel_heights
first_col = cowplot::plot_grid(p2, p1, ncol = 1, rel_heights = c(1, 3),align = "v",axis = "l")
first_col_bar <- cowplot::plot_grid(p2, p1bar, ncol = 1, rel_heights = c(1, 3),align = "v",axis = "l")
if (numpairs>=1) {
  first_col_n = cowplot::plot_grid(p2n, p1, ncol = 1, rel_heights = c(1, 3),align = "v",axis = "l")
  first_col_nuncor <- cowplot::plot_grid(p2nuncor, p1, ncol = 1, rel_heights = c(1, 3),align = "v",axis = "l")
}
first_col_pseudolog <- cowplot::plot_grid(p2pseudolog, p1pseudolog, ncol = 1, rel_heights = c(1, 3),align = "v",axis = "l")
first_col_pseudolog_bar <- cowplot::plot_grid(p2pseudolog, p1pseudolog_bar, ncol = 1, rel_heights = c(1, 3),align = "v",axis = "l")
first_col_bar
first_col_pseudolog_bar
  
Warning <- cowplot::ggdraw() + cowplot::draw_label("Warning!\nAxes not \nat same \nscale!", fontface='bold')
warninglog <- cowplot::ggdraw() + cowplot::draw_label("X axis \nlog 10 \ntransformed!", fontface='bold')
warningpseudolog <- cowplot::ggdraw() + cowplot::draw_label("X axis \npseudo-log 10 \ntransformed!", fontface='bold')
second_col = cowplot::plot_grid(legend, p3, ncol = 1, rel_heights = c(1, 3))
if (numpairs>=1) {
  second_col_n = cowplot::plot_grid(Warning, p3n, ncol = 1, rel_heights = c(1, 3))
  second_col_nuncor = cowplot::plot_grid(Warning, p3nuncor, ncol = 1, rel_heights = c(1, 3))
}
second_col_log = cowplot::plot_grid(legend, p3, ncol = 1, rel_heights = c(1, 3))
second_col_pseudolog = cowplot::plot_grid(legend, p3, ncol = 1, rel_heights = c(1, 3))

perfect = cowplot::plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1.2),axis="l")
if (numpairs>=1) {
  perfect_n = cowplot::plot_grid(first_col_n, second_col_n, ncol = 2, rel_widths = c(3, 1.2),axis="l")
  perfect_nuncor = cowplot::plot_grid(first_col_nuncor, second_col_nuncor, ncol = 2, rel_widths = c(3, 1.2),axis="l")
}
perfect_pseudolog <- cowplot::plot_grid(first_col_pseudolog, second_col_pseudolog, ncol = 2, rel_widths = c(3, 1.2),axis="l")
perfect_bar <- cowplot::plot_grid(first_col_bar, second_col, ncol = 2, rel_widths = c(3, 1.2),axis="l")
perfect_bar_pseudolog <- cowplot::plot_grid(first_col_pseudolog_bar, second_col_pseudolog, ncol = 2, rel_widths = c(3, 1.2),axis="l")
perfect
if (numpairs>=1) {
  perfect_n
  perfect_nuncor
}
perfect_pseudolog
perfect_bar
perfect_bar_pseudolog
  
if (is.null(ylims.plot)) {
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,".pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect); dev.off()
    if (numpairs>=1) {
      pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-n.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_n); dev.off()
      pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-nuncor.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_nuncor); dev.off()
    }
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-xpseudolog10.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_pseudolog); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-bar.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_bar); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-bar-xpseudolog10.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_bar_pseudolog); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-bar-small.pdf",sep=""),sep = '')),height=8,width=2,useDingbats=F); print(first_col_bar); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-bar-xpseudolog10-small.pdf",sep=""),sep = '')),height=8,width=2,useDingbats=F); print(first_col_pseudolog_bar); dev.off()
} else if (!is.null(ylims.plot)) {
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],".pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect); dev.off()
    if (numpairs>=1) {
      pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-n.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_n); dev.off()
      pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-nuncor.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_nuncor); dev.off()
    }
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-xpseudolog10.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_pseudolog); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-bar.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_bar); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-bar-xpseudolog10.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(perfect_bar_pseudolog); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-bar-small.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(first_col_bar); dev.off()
    pdf(file.path(PathToParameterPlots,paste(paste0(filename,"-y",ylims.plot[1],"-",ylims.plot[2],"-bar-xpseudolog10-small.pdf",sep=""),sep = '')),height=8,width=8,useDingbats=F); print(first_col_pseudolog_bar); dev.off()
}

closeAllConnections()
dev.off()
rm(list=ls())
gc()

