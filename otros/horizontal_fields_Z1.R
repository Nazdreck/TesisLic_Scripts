########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ncdf4)          # For manipulate .nc files
require(ggplot2)        # Data Visualisation
require(metR)           # Tools for Meteorological Analysis
require(colorspace)     # Paletas de Colores.
require(abind)          # Combine Multidimensional Arrays
require(ggpubr)         # "ggplot2" Based Publication Ready Plots
require(maps)           # Draw geographical maps
require(mapproj)        # Map Projections
require(ggtext)
require(gtools)

########## Predefiniendo Variables ############################################
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(0,-85,-2.5); ny = length(lats)
plvl = 300
days = 1:92; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "1900-06-01")), 6,10)
stats = c("2002", "2019", "climatology"); ns = length(stats)

DATA = array(
  double(nx*ny*nt*ns), c(nx,ny,nt,ns), 
  dimnames = list(lons = lons, lats = lats, dates = dates, stats = stats)
)

########## Extrayendo datos ###################################################

for (s in stats) {
  nc = nc_open(
    here("inputs/era5_plvls", paste0("gh_",s,".nc"))
  )
  data_extracted = ncvar_get(
    nc,
    start = c(
      1,
      which(ncvar_get(nc, "lats") == lats[1]),
      which(ncvar_get(nc, "plvl") == plvl),
      which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "07-01")
    ),
    count = c(
      -1,ny,1,nt
    )
  )
  
  DATA[,,,s] = data_extracted
  nc_close(nc); rm(nc)
}
rm(data_extracted)

########## Analisis Espectral #################################################
FOURIER = apply(DATA/9.8, 2:4, FilterWave, k = 1)
FOURIER = apply(FOURIER, c(1,2,4), mean)

########## Acomodando datos ###################################################
ANOM = FOURIER[,,c("2002","2019")]
for (i in 1:2) {
  ANOM[,,i] = ANOM[,,i] - FOURIER[,,"climatology"]
}

data = abind(ANOM, FOURIER[,,"climatology"], along = 3)
########## Elementos para graficar ############################################
brks_clim = seq(-140,140,20); nbc = length(brks_clim)
brks_anom = seq(-60,60,10); nba = length(brks_anom)

cols_clim = c(
  sequential_hcl((nbc-1)/2, palette = "Blues", rev = FALSE),
  sequential_hcl((nbc-1)/2, palette = "Reds", rev = TRUE)
)

cols_anom = c(
  sequential_hcl((nba-1)/2, palette = "Blues", rev = FALSE),
  sequential_hcl((nba-1)/2, palette = "Reds", rev = TRUE)
)

brks = list(brks_anom,brks_anom,brks_clim)
cols = list(cols_anom,cols_anom,cols_clim)

ttl = list(
  expression("Anomalía de Z"[1]~"[m] (2002)"),
  expression("Anomalía de Z"[1]~"[m] (2019)"),
  expression("Climatología de Z"[1]~"[m] (1979-2018, s/2002)")
)

# ttl = list("(a)", "(b)", "(c)")

##### theming
# x axis
tm.axis_text_x_bottom = list("top" = element_text(), "center" = element_blank(), "bot" = element_blank())
tm.axis_text_x_top = list("top" = element_blank(), "center" = element_blank(), "bot" = element_text())
tm.axis_ticks_x_bottom = list("top" = element_line(), "center" = element_line(), "bot" = element_blank())
tm.axis_ticks_x_top = list("top" = element_blank(), "center" = element_line(), "bot" = element_line())

########## Graficos ###

GG = vector("list", 3); names(GG) = stats
for (s in 1:ns) {
  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    geop = array(data[,,s])
  )
  
  GG[[s]] = ggplot(df) +
    aes(x = lons, y = lats) +
    geom_contour_fill(aes(z = geop), breaks = brks[[s]]) +
    scale_fill_gradientn(
      breaks = brks[[s]],
      colors = cols[[s]],
      limits = range(brks[[s]])
    ) +
    guides(fill = "none") +
    geom_contour2(aes(z = geop), breaks = brks[[s]], color = "gray70", alpha = 0.6, size = 0.3) +
    geom_text_contour(aes(z = geop), breaks = brks[[s]], stroke = 0.15, size = 2.5) +
    geom_map(data = map_data("world2"), aes(x = long, y = lat, map_id = region), map = map_data("world2"), fill = NA, color = "black", size = 0.1) +
    scale_x_longitude(breaks = seq(60,300,60), limits = range(lons)) +
    scale_y_latitude(breaks = seq(-75,-15,15), limits = c(-85,-5)) +
    labs(title = ttl[[s]]) +
    theme_light() +
    theme(
      plot.margin = margin(0,0,0,0, unit = "cm"),
      plot.title = element_text(size = 6, vjust = 0, hjust = 0, face = "bold", margin = margin(0,0,0,0, unit = "cm")),
      axis.title = element_blank(),
      axis.text = element_text(size = 5)
    )
}

figure = ggarrange(
  plotlist = GG, ncol = 1, nrow = 3
)

ggsave(
  here("outputs/ftropos", "Z_fields.png"),
  figure, device = "png", width = 12, height = 9, units = "cm", dpi = 200
)























