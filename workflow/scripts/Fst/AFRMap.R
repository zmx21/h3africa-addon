library(sf)
library(raster)
library(dplyr)
library(spData)
library(tmap)
library(leaflet)
library(cartogram)
library(ggplot2)
library(ggrepel)
library(grid)
library(ggplotify)
library(gridExtra)
library(ggnetwork)
library(latex2exp)
library(grConvert)
out_path <- '../results/Fst/'


#Download latitude and longtitude of 1KG Populations
#https://www.internationalgenome.org/data-portal/population
igsr_locations <- data.table::fread(glue::glue('{out_path}/igsr_populations.tsv')) %>% dplyr::filter(`Population code` %in% c("YRI","LWK","GWD","MSL","ESN")
) %>% dplyr::select(POP = `Population code`,Longtitude = `Population longitude`,Latitude = `Population latitude`) %>% dplyr::mutate(Color = '1KG/AFGR')
locations <- rbind(data.frame(POP = 'TB-DAR',Longtitude = 39.279556,Latitude = -6.802353,Color = 'TB-DAR'),igsr_locations)
AFGR_Locations <- data.frame(POP = c('Uganda','Ethiopia','Namibia','Zulu'),Longtitude = c(32.55,38.7,17.083333,30.383333),Latitude = c(0.316666667,9.033333333,-22.56666667,-29.616667),Color = 'AFGR Only')
locations <- rbind(locations,AFGR_Locations)
locations$POP <- as.character(locations$POP)
africa_pop = world %>% 
  filter(continent %in%  c("Africa"))

#Construct Comb Matrix
fst_results <- data.table::fread(glue::glue("{out_path}TBDAR.1KG.maf0p05.details.Fst"),skip = 1,header = F) %>% dplyr::filter(V2 != 'ACB' & V3 != 'ACB')
AFR_POP <- unique(fst_results$V2)
comb_matrix <- matrix(NA,nrow = length(AFR_POP),ncol = length(AFR_POP))
rownames(comb_matrix) <- AFR_POP
colnames(comb_matrix) <- AFR_POP

for(i in 1:nrow(fst_results)){
  cur_row <- which(rownames(comb_matrix) == fst_results$V2[i])
  cur_col <- which(colnames(comb_matrix) == fst_results$V3[i])
  comb_matrix[cur_row,cur_col] <- fst_results$V4[i]
  comb_matrix[cur_col,cur_row] <- fst_results$V4[i]
  
}
rownames(comb_matrix)[rownames(comb_matrix) == 'TBDAR'] <- 'TB-DAR'
colnames(comb_matrix)[colnames(comb_matrix) == 'TBDAR'] <- 'TB-DAR'


non_tb_dar_loc <- dplyr::filter(locations,POP != 'TB-DAR' & Color != 'AFGR Only') 
edges <-non_tb_dar_loc %>% dplyr::select(xend = Longtitude,yend =  Latitude,POP=POP)
edges$Fst <- comb_matrix[non_tb_dar_loc$POP,'TB-DAR']
edges$x <- locations$Longtitude[locations$POP == 'TB-DAR']
edges$y <- locations$Latitude[locations$POP == 'TB-DAR']


p1 <- ggplot(africa_pop) + 
  geom_sf(aes(geometry = geom)) + 
  geom_point(aes(x = Longtitude, y = Latitude),data = locations) + 
  geom_label_repel(data = locations %>% dplyr::filter(!POP%in%c('YRI','TB-DAR','MSL') & Color != 'AFGR Only'),mapping = aes(label = as.character(POP),x = Longtitude, y = Latitude,fill = factor(Color)),size = 3,nudge_x = 8,nudge_y = 5,direction = 'both') +
  geom_label_repel(data = locations %>% dplyr::filter(POP=='YRI'),mapping = aes(label = as.character(POP),x = Longtitude, y = Latitude,fill = factor(Color)),size = 3,nudge_x = 0,nudge_y = 5,direction = 'both') +
  geom_label_repel(data = locations %>% dplyr::filter(POP=='TB-DAR'),mapping = aes(label = as.character(POP),x = Longtitude, y = Latitude,fill = factor(Color)),size = 3,nudge_x = 0,nudge_y = -5,direction = 'both') + 
  geom_label_repel(data = locations %>% dplyr::filter(POP=='MSL'),mapping = aes(label = as.character(POP),x = Longtitude, y = Latitude,fill = factor(Color)),size = 3,nudge_x = -2,nudge_y = -5,direction = 'both') + 
  geom_label_repel(data = locations %>% dplyr::filter(!Color %in% c('TB-DAR','1KG/AFGR')),mapping = aes(label = as.character(POP),x = Longtitude, y = Latitude,fill = factor(Color)),size = 3,nudge_y = 5,nudge_x = 1,direction = 'both') +
  geom_edges(aes(x=x,y=y,xend=xend,yend = yend,color = Fst),data = edges[c(1,2,5),],curvature = 0.2,size=1.1)+
  geom_edges(aes(x=x,y=y,xend=xend,yend = yend,color = Fst),data = edges[c(3,4),],curvature = -0.2,size=1.1)+ 
  scale_color_gradient(low="blue", high="red") + theme(plot.margin = margin(t = 1, r = 0, b = 0, l = 0.5,'cm')) + labs(color = TeX("$F_{ST}$")) + guides(
    fill = guide_legend(
      title = "Source",
      override.aes = aes(label = "")
    )
  )



comb_matrix_merged_reorder <- comb_matrix[c('GWD','MSL','YRI','ESN','TB-DAR','LWK'),c('GWD','MSL','YRI','ESN','TB-DAR','LWK')]
col2 <- colorRampPalette(c("#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
pdf(file = '../results/Plots/corrplot.pdf',width = 5,height = 5)
corrplot::corrplot(signif(comb_matrix_merged_reorder,2),
                   method = 'shade',type = 'upper',
                   is.corr = F,number.digits=4,diag = T,cl.pos = 'n', cl.lim = c(0,0.015),col = col2(200),
                   addCoef.col='black',tl.col = 'black',na.label = NA)
dev.off()

grConvert::convertPicture('../results/Plots/corrplot.pdf', '../results/Plots/corrplot.svg')
corr_plot <- grImport2::pictureGrob(grImport2::readPicture('../results/Plots/corrplot.svg'))

cowplot::plot_grid(p1,corr_plot,ncol = 2,labels = c('A)','B)'))
