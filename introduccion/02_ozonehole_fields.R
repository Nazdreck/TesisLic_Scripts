########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ncdf4)          # For manipulate .nc files
require(ggplot2)        # Data Visualisation
require(metR)           # Tools for Meteorological Analysis
require(colorRamps)     # Builds color Tables
require(RColorBrewer)   # Paletas de Colores.
require(colorspace)     # Paletas de Colores.
require(abind)          # Combine Multidimensional Arrays
require(ggpubr)         # "ggplot2" Based Publication Ready Plots
require(maps)           # Draw geographical maps
require(mapproj)        # Map Projections

source(here("functions/gen_breaks.R"))

########## Extraccion de datos ################################################
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(90,-90,-2.5); ny = length(lats)
years = c(2001:2003, 2018:2020); na = length(years)
days = c(27, 19); nd = length(days)
days_rep = rep(days, each = 3); ndr = length(days_rep)

dates = paste(days_rep,"SEP",years, sep = "-")

DATA = NULL

nc = nc_open(
  here("inputs/era5_ozone/tco_09-10_79-19.nc")
)

for (a in 1:(na-1)) {
  if (a <= 3) d = days[1]
  if (a > 3) d = days[2]
  data_extracted = ncvar_get(
    nc,
    varid = "tco",
    start = c(
      1,1,d,
      which(ncvar_get(nc, "years") == years[a])
    ),
    count = c(-1,-1,1,1)
  )
  DATA = abind(DATA,data_extracted, along = 3)
}

nc_close(nc);rm(nc)

nc = nc_open(
  here("inputs/era5_ozone/tco_09_2020_day.nc")
)

data_extracted = ncvar_get(
  nc,
  varid = names(nc$var)[2],
  start = c(1,1,days[2]),
  count = c(-1,-1,1)
) 

DATA = abind(DATA, data_extracted, along = 3)

nc_close(nc);rm(nc);rm(data_extracted)

dimnames(DATA) = list(lons = lons, lats = lats, dates = dates)
########## Preparando los datos ###############################################
dobson = 2.1415e-5

d360 = DATA["180",,]
d360 = abind(d360,DATA[which(lons > 180),,], along = 1)
d360 = abind(d360, DATA[which(lons <= 180),,], along = 1)
lons = seq(-180,180,2.5); nx = length(lons)
dimnames(d360) = list(lons = lons, lats = lats, dates = dates)

d180 = d360[,which(lats <= 0),]/dobson
lats = seq(0,-90,-2.5); ny = length(lats)

########## Elementos para graficar ############################################

brks = seq(100,480,5)
ref = seq(100,480,40)
labs = rep("", length(brks))
for (i in 1:length(brks)) {
  if (any(as.character(brks[i]) == as.character(ref))) {
    labs[i] = brks[i]
  }
}

########## Figura #############################################################
GG = vector("list", 6)

for (a in 1:na) {
  df = data.frame(
    lons = rep(lons,ny),
    lats = rep(lats, each = nx),
    ozone = array(d180[,,a], nx*ny)
  )
  
  print(paste0("Graficando ", dates[a], " ..."))
  
  GG[[a]] = ggplot(df) +
    aes(x = lons, y = lats, z = ozone) +
    geom_contour_fill(breaks = brks, show.legend = TRUE) +
    scale_fill_gradientn(
      colors = sequential_hcl(length(brks)-1, palette = "inferno"), 
      limits = range(brks), 
      breaks = brks, 
      labels = labs 
    ) +
    guides(
      fill = guide_colourstrip(
        title = "UD", barwidth = 14, barheight = 0.3, title.position = "right", title.hjust = 15,
        label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
        )
    ) +
    geom_map(data = map_data("world"), map = map_data("world"), aes(map_id = region), fill = NA, color = "white", size = 0.1, inherit.aes = FALSE) +
    coord_map("ortho", orientation = c(-90,0,0)) +
    geom_contour2(breaks = c(220), linetype = "dashed", size = 0.2, color = "deepskyblue") + 
    scale_y_continuous(breaks = c(-60,-30), labels = NULL, name = NULL) +
    scale_x_continuous(breaks = seq(-120,180,60), labels = NULL, name = NULL) +
    annotate("text", x = -120, y = -30, label = "30°S", color = "white", size = 1.75) +
    annotate("text", x = -120, y = -60, label = "60°S", color = "white", size = 1.75) +
    labs(title = dates[a]) +
    theme(
      panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.05, linetype = "dashed", color = "gray15"),
      panel.background = element_blank(), axis.ticks=element_blank(),
      legend.position="bottom",plot.title = element_text(hjust = -0.1, vjust = -5, size = 7),
      legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold")
      )
  
  print("Hecho.")
}

print("Combinando figuras...")
panel = ggarrange(plotlist = GG, nrow = 2, ncol = 3, common.legend = TRUE, legend = "bottom")
print("Hecho.")

print("Exportando figuras...")
ggsave(
  here("outputs/introduccion", "ozone_hole2.png"),
  panel, device = "png", width = 16, height = 12, units = "cm", dpi = 200
)
print("Hecho.")



