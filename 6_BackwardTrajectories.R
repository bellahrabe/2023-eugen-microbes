basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
package.list <- setdiff(package.list,basic.packages)
if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
rm(list=ls()); closeAllConnections(); gc()

# # this reads in the output of hysplit, creates shp, csv and rds files per cruise and writes them out
# # also creates plots of x-day back trajectories per cruise and for all cruises together

if (!require('ggplot2')) {install.packages('ggplot2')}; library(ggplot2)
if (!require('dplyr')) {install.packages('dplyr')}; library(dplyr)


folder.names <- list.dirs(path="C:/Users/i.hrabedeangelis/Dropbox/B/sci/phd/v01/bt_raw",full.names = T)
folder.names <- folder.names[-c(1)]
folder.names.short <- list.dirs(path="C:/Users/i.hrabedeangelis/Dropbox/B/sci/phd/v01/bt_raw",full.names = F) # including BT start height
folder.names.short <- folder.names.short[-c(1)]
folder.names.very.short <- sub('.*_', '', folder.names.short) # cruise name only


folderforresults <- file.path("C:/Users/i.hrabedeangelis/Dropbox/B/sci/phd/v01/bt_out_new")
dir.create(folderforresults)
dir.create(paste0(folderforresults,"/bycruise/"))

Ilikethismanydays <- 3 # only for plots, not for data output, length of BTs wanted
RainCutoff <- 0.5 # mm per h
#BTparamOfInterest <- "timedif" # choose between timedif, heightaboveground, pressure, theta, air_temp, rainfall, mixdepth, relhumid, H2Omixra,sunflux
BTparamsOfInterest <- c("timedif","heightaboveground","pressure","theta","air_temp","rainfall","mixdepth","relhumid","H2Omixra","sunflux", paste0("has",RainCutoff,"mm_p_h_rainfall_withinPast",Ilikethismanydays,"days"))

# #  plot prep #################################################################################
my_colors <- colorspace::sequential_hcl(n=7,palette="Plasma")
theme_set(theme_bw())

countries <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
coastlines <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")

# data.all <- matrix(ncol=22,nrow=0)
# colnames(data.all) <- c("NA","NA2","year","month","day","hour","minute","sec","timedif","lat","lon","heightaboveground","pressure","theta","air_temp","rainfall","mixdepth","relhumid","H2Omixra","sunflux","index","cruise")


# # goes through all folders (one per cruise), reads in the back trajectories, creates a shape file per cruise and writes it out as shape file and rds
for (k in 1:length(folder.names)) { 
  file.names <- list.files(path = folder.names[k] ,full.names = T)
  data.main <- matrix(ncol=22,nrow=0)
  colnames(data.main) <- c("NA","NA2","year","month","day","hour","minute","sec","timedif","lat","lon","heightaboveground","pressure","theta","air_temp","rainfall","mixdepth","relhumid","H2Omixra","sunflux","index","cruise")
  
  # # read in the data for the current cruise
  for (i in 1:length(file.names)) {
    file <- file.names[i]
    lines <- readLines(file)
    # Identify the indices of potential data rows with 20 numeric values
    potential_data_indices <- sapply(strsplit(lines, "\\s+"), function(x) sum(!is.na(as.numeric(x))) == 20)
    # Find the start index of the trajectory data by locating the first TRUE value
    start_index <- which(potential_data_indices)[1]
    
    data_lines <- lines[start_index:length(lines)]  # Read data
    data <- read.table(text=data_lines, header = FALSE, stringsAsFactors = FALSE)
    
    data$index <- i
    data$cruise <- folder.names.very.short[k]
    # Create a data frame for each trajectory
    colnames(data) <- c("NA","NA2","year","month","day","hour","minute","sec","timedif","lat","lon","heightaboveground","pressure","theta","air_temp","rainfall","mixdepth","relhumid","H2Omixra","sunflux","index","cruise")
    
    data.main <- rbind(data.main,data)
  }
  rm(data)
  
  # data.all <- rbind(data.all,data.main)
  # Step 4: write cruisefile as rds and csv
  saveRDS(data.main, paste0(folderforresults,"/bycruise/",folder.names.short[k],".rds",sep=""))
  write.csv(data.main, paste0(folderforresults,"/bycruise/",folder.names.short[k],".csv",sep=""), row.names = FALSE)
  
  # # cut off after certain time
  data.main <- data.main %>% filter(timedif>=-(Ilikethismanydays*24)) ## 
  
  
  # # cut bts if rain occured  and adds two columns (1 whether it was cut due to rainfall, 1 at what time (in h)#####
  data.main <- data.main %>%
    group_by(index) %>%
    mutate(has_high_rainfall = any(rainfall >= RainCutoff),
           lastrainfall = ifelse(has_high_rainfall, max(timedif[rainfall >= RainCutoff]), NA)) %>%
    ungroup() %>%
    filter(is.na(lastrainfall) | timedif > lastrainfall) # %>%
  #select(-has_high_rainfall, -max_timedif)
  colnames(data.main)[ncol(data.main)-1] <-  paste0("has",RainCutoff,"mm_p_h_rainfall_withinPast",Ilikethismanydays,"days")
  
  
  # Step 3: Create a Geometry Column
  sf_data <- sf::st_as_sf(data.main, coords = c("lon", "lat"), crs = 4326)
  # Step 4: Create Shapefile
  # st_write(sf_data, paste0(folderforresults,"/bycruise/",folder.names.short[k],".shp",sep=""))
  
  # # create sf line object to check whether lines cross with land masses
  subset.bts.allsample.lines.sf<-sf_data %>% group_by(index) %>% summarize(do_union=FALSE) %>% sf::st_cast("LINESTRING")
  # # Check for intersections between line elements and coastlines. Reurns which element in each BT interacts with land
  line_coastline_intersection <- sf::st_intersects(subset.bts.allsample.lines.sf, coastlines)
  # Calculate the number of lines crossing coastlines
  lines_crossing_count <- sum(sapply(line_coastline_intersection, any))
  # Calculate the percentage of lines crossing coastlines
  lines_crossing_percentage <- (lines_crossing_count / nrow(subset.bts.allsample.lines.sf))
  lines_crossing_TF <- lines_crossing_percentage>0
  
  for (i in 1:length(BTparamsOfInterest)) {
    cur_param <- BTparamsOfInterest[i]
    # # PLOTS normal
    pdf(file.path(paste(folderforresults,"/bycruise/",folder.names.short[k],"_",cur_param,"_",RainCutoff,"raincut_",Ilikethismanydays,"day_bt.pdf",sep="")),height=7.5,width=7.5,useDingbats=F)
    p2 <- ggplot(sf_data) +
      geom_sf(data=countries, fill="grey",color="grey",linewidth=0.01)+
      geom_sf(data=coastlines, colour="black",linewidth=0.1)+
      geom_sf(aes(color = !!sym(BTparamsOfInterest[i])),size=0.5,shape=16)  ## parameters from bt: heightaboveground, pressur, theta, air_tmp, rainfll, mixdpth, relhumd, H2=mixr, sunflux
    if (is.numeric(sf_data[[cur_param]])) {
      p2 <- p2 +
        scale_color_gradientn(colors=my_colors) 
    } else {
      p2 <- p2 +
        guides(color = guide_legend(override.aes = list(size = 2)))  
    }
    p2 <- p2 +
      coord_sf(xlim = c(-55,35), ylim = c(-10, 85), expand = F) +
      ggtitle(paste0(folder.names.short[k]," - raincutoff ",RainCutoff," mm per h - ",Ilikethismanydays,"-day back trajectories, 200 m",sep=""))
    print(p2)
    dev.off()
    
    # # PLOTS zoom
    data_extent <- sf::st_bbox(sf_data)
    pdf(file.path(paste(folderforresults,"/bycruise/",folder.names.short[k],"_",cur_param,"_",RainCutoff,"raincut_",Ilikethismanydays,"day_bt_zoom.pdf",sep="")),height=7.5,width=7.5,useDingbats=F)
    p2 <- ggplot(sf_data) +
      geom_sf(data=countries, fill="grey",color="grey",linewidth=0.01)+
      geom_sf(data=coastlines, colour="black",linewidth=0.1)+
      geom_sf(aes(color = !!sym(BTparamsOfInterest[i])),size=0.5,shape=16)  ## parameters from bt: heightaboveground, pressur, theta, air_tmp, rainfll, mixdpth, relhumd, H2=mixr, sunflux
    if (is.numeric(sf_data[[cur_param]])) {
      p2 <- p2 +
        scale_color_gradientn(colors=my_colors) 
    } else {
      p2 <- p2 +
        guides(color = guide_legend(override.aes = list(size = 2)))  
    }
    p2 <- p2 +
      coord_sf(xlim = c((data_extent$xmin)-5, (data_extent$xmax)+5),
               ylim = c((data_extent$ymin)-5, (data_extent$ymax)+5),
               expand = FALSE) +
      ggtitle(paste0(folder.names.short[k]," - raincutoff ",RainCutoff," mm per h - ",Ilikethismanydays,"-day back trajectories, 200 m",sep=""))
    print(p2)
    dev.off()
  } 
}
