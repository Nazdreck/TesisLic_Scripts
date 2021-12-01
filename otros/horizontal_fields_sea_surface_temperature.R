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

########## Extrayendo datos ###################################################
### Nino 3.4
NINO34 = read.table(
  here("inputs/ENSO/nino34.txt"), sep = "", dec = ".", row.names = 1
)

NINO34 = array(t(as.matrix(NINO34)))


### SST
lons = seq(0,358, 2); nx = length(lons)
lats = seq(88,-88,-2); ny = length(lats)
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = as.character(1981:2019); na = length(years)

nc = nc_open(here("inputs/ENSO/sst.mnmean.nc"))

# months = c(
#   which(as.character(as.Date(ncvar_get(nc, "time"), origin = "1800-01-01")) == "1981-01-01"):
#   which(as.character(as.Date(ncvar_get(nc, "time"), origin = "1800-01-01")) == "2019-12-01")
# ); nt = length(months)
# dates = as.character(as.Date(ncvar_get(nc, "time")[months], origin = "1800-01-01"))

SST = ncvar_get(
  nc,
  varid = names(nc$var)[2],
  start = c(
    1,
    1,
    which(as.character(as.Date(ncvar_get(nc, "time"), origin = "1800-01-01")) == "1981-01-01")
  ),
  count = c(
    -1, -1, nm*na
  )
)
SST = array(array(SST), c(nx,ny,nm,na), list(lons = lons, lats = lats, months = months, years = years))
########## Preparando SST #####################################################
# # Quito tendencia lineal
# RL = array(
#   rep(NA,nx*ny*2), c(nx,ny,2),
#   dimnames = list(lons = lons, lats = lats, param = c("pend", "ordorg"))
# )
# 
# ALL = array(vector("logical", nx*ny), c(nx,ny))
# ANY = array(vector("logical", nx*ny), c(nx,ny))
# 
# for (y in 1:ny) {
#   for (x in 1:nx) {
#     ALL[x,y] = all(is.na(SST[x,y,]))
#     ANY[x,y] = any(is.na(SST[x,y,]))
#   }
# }
# 
# 
# for (y in 1:ny) {
#   for (x in 1:nx) {
#     if (!ANY[x,y]) {
#       RL[x,y,"pend"] = summary(lm(SST[x,y,] ~ months))$coef[2]
#       RL[x,y,"ordorg"] = summary(lm(SST[x,y,] ~ months))$coef[1]
#     }
#   }
# }
# 
# SST_RL = array(
#   double(nx*ny*nt), c(nx,ny,nt),
#   dimnames = list(lons = lons, lats = lats, dates = dates)
# )
# for (t in 1:nt) {
#   for (y in 1:ny) {
#     for (x in 1:nx) {
#       SST_RL[x,y,t] = RL[x,y,"pend"]*months[t]+RL[x,y,"ordorg"]
#     }
#   }
# }
# 
# # CLimatologia
# CLIM = apply(SST_RL[,,1:360], c(1,2), mean)
# 
# 
# # Extraigo meses de interes
# tt = c(-1,0,1)
# mm = c(
#   dates[which(dates == "2002-08-01") + tt],
#   dates[which(dates == "2019-08-01") + tt]
# )
# ntt = length(mm)
# DATA = SST_RL[,,mm]
# 
# 
# ### Anomalia
# ANOM = array(double(nx*ny*ntt), c(nx,ny,ntt), dimnames = list(lons = lons, lats = lats, months = mm))
# for (t in 1:ntt) {
#   ANOM[,,t] = DATA[,,t] - CLIM 
# }

#########################
CLIM = apply(SST, c(1,2,3), mean)

mm = c("07", "08", "09")
yy = c("2002", "2019")
DATA = SST[,,mm, yy]

ANOM = array(double(nx*ny*3*2), c(nx,ny,3,2), dimnames = list(lons = lons, lats = lats, months = mm, years = yy))
for (a in 1:2){
  for (m in mm) {
    ANOM[,,m,a] = DATA[,,m,a] - CLIM[,,m]
  }
}
########## Elementos para graficar ############################################
brks = seq(-5,5,0.5); nb = length(brks)
labs = as.character(brks)
labs[which(odd(1:nb))] = " "

cols = c(
  sequential_hcl((nb-1)/2, palette = "Blues", rev = FALSE),
  sequential_hcl((nb-1)/2, palette = "YlOrRd", rev = TRUE)
)

A = matrix(
  c(rep(yy,each = 3), rep(mm, 2)), ncol = 2, nrow = 6
)

########## Graficos ###########################################################
########## SST
GG = vector("list", 6)

brks = c(-5,seq(-3,3,0.5),5); nb = length(brks)
labs = as.character(brks)
labs[which(odd(1:nb))] = " "

A = cut(array(ANOM[,,m,a]), breaks = brks, labels = labs)

for (i in 1:6) {
  a = A[i,1]
  m = A[i,2]
  
  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    sst = array(ANOM[,,m,a])
  )
  
  GG[[i]] = ggplot(df) +
    aes(x = lons, y = lats, fill = sst) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradientn(
      colors = cols, na.value = NA,
      breaks = brks,
      labels = labs,
      limits = range(brks)
    ) +
    guides(
      fill = guide_colorstrip(
        title = "SST [K]  ", ticks = TRUE, ticks.color = "gray30", barwidth = 10, barheight = 0.3, label.vjust = 2
      )
    ) +
    geom_map(data = map_data("world2"), aes(x = long, y = lat, map_id = region), map = map_data("world2"), fill = "gray60", color = "black", size = 0.05) +
    #coord_map(projection = "mercator", orientation = c(-90,0,0), ylim = c(-90,0)) +
    scale_x_longitude(breaks = seq(0,360,60)) +
    scale_y_latitude(breaks = seq(-60,60,30), limits = c(-80,80)) +
    labs(title = paste0(m,"-",a)) +
    theme_light() +
    
    theme(
      plot.title = element_text(size = 6, face = "bold", vjust = -3),
      plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0, unit = "cm"),
      legend.title = element_text(size = 5, face = "bold"),
      legend.text = element_text(size = 4),
      legend.margin = margin(c(0,0,0,0), unit = "cm"),
      legend.box.spacing = unit(0, units = "cm"),
      axis.title = element_blank(),
      axis.text = element_text(size = 5)
    ) +
    coord_fixed(1) 
}

figure = ggarrange(plotlist = GG, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/ftropos/sst_fields.png"),
  figure, device = "png", width = 16, height = 6.75, units = "cm", dpi = 300
)

########## NINO3.4
brks = seq(-2.6,2.6,0.1); nb = length(brks)
# labs = as.character(brks)
# labs[which(odd(1:nb))] = " "

cols = c(
  sequential_hcl((nb-1)/2, palette = "Blues", rev = FALSE),
  sequential_hcl((nb-1)/2, palette = "YlOrRd", rev = TRUE)
)

df = data.frame(
  time = 1:360,
  var = NINO34
)

GG = ggplot(df) +
  aes(x = time, y = var) +
  geom_line(size = 0.25, alpha = 0.7) +
  scale_y_continuous(name = "Anomalía [K]", breaks = -3:3, limits = c(-3,3), expand = c(0,0)) +
  scale_x_continuous(name = NULL, breaks = c(seq(1,360,12*5),360), labels = seq(1990,2020,5), limits = range(1:360), expand = c(0,0)) +
  geom_hline(yintercept = 1, color = "red", size = 0.25, alpha = 0.7) +
  geom_hline(yintercept = -1, color = "blue", size = 0.25, alpha = 0.7) +
  geom_vline(xintercept = 154, size = 0.25, alpha = 0.7) +
  annotate("text", x = 11, y = 1.2, label = "El Niño", color = "red", size = 1.5) +
  annotate("text", x = 11, y = -1.2, label = "La Niña", color = "blue", size = 1.5) +
  labs(title = "El Niño 3.4") +
  theme_light() +
  theme(
    plot.margin = margin(t = 0.1, r = 0.3, b = 0.1, l = 0.1, unit = "cm"),
    plot.title = element_text(face = "bold", size = 6, vjust = -3),
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 6)
  )


ggsave(
  here("outputs/ftropos", "nino34.png"),
  GG, device = "png", width = 12, height = 5.0625, units = "cm", dpi = 200
)






















