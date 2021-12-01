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
require(gtools)


source(here("functions/gen_breaks.R"))

########## Predefiniendo variables ############################################
vars = c("gh", "pv"); nv = length(vars)
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats_hs = seq(-15, -90, -2.5); lats_hn = seq(90,15,-2.5); ny = length(lats_hs)
days = 1:31; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)

SH = array(
  double(nx*ny*nv), dim = c(nx,ny,nv), 
  dimnames = list(lons = lons, lats = lats_hs, vars = vars)
)

NH = array(
  double(nx*ny*nv), dim = c(nx,ny,nv), 
  dimnames = list(lons = lons, lats = lats_hn, vars = vars)
)

for (v in 1:nv){
  nc = nc_open(
    here("inputs/era5_plvls", paste0(vars[v], "_climatology.nc"))
  )
  SH[,,v] = apply(
    ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == lats_hs[1]),
        which(ncvar_get(nc, "plvl") == 10),
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "08-01")
      ), 
      count = c(-1,-1,1,nt)
    ),
    c(1,2), mean
  )
  
  NH[,,v] = apply(
    ncvar_get(
      nc,
      start = c(
        1,
        1,
        which(ncvar_get(nc, "plvl") == 10),
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "01-01")
      ), 
      count = c(-1,ny,1,nt)
    ),
    c(1,2), mean
  )
  nc_close(nc); rm(nc)
}

##########  Preparando datos ##################################################
#####  Latitudes
lats = list(SH = lats_hs,NH = lats_hn)
#####  Datos
DATA = list(SH = SH, NH = NH)

for (i in 1:2) {
  lons = seq(0, 357.5, 2.5); nx = length(lons)
  d360 = DATA[[i]]["180",,]
  d360 = abind(d360, DATA[[i]][which(lons > 180),,], along = 1)
  d360 = abind(d360, DATA[[i]][which(lons <= 180),,], along = 1)
  lons = seq(-180,180,2.5); nx = length(lons)
  dimnames(d360) = list(lons = lons, lats = lats[[i]], vars = vars)
  DATA[[i]] = d360
}


########## Elementos para graficar ############################################
brks = seq(-700,700,50); nb = length(brks)
labs = as.character(brks)
labs[which(even(1:nb))] = " "
gbrks = seq(26, 31, 0.5)

pv_cols = c(
  sequential_hcl((nb-1)/2, palette = "GnBu", rev = FALSE),
  sequential_hcl((nb-1)/2, palette = "YlOrRd", rev = TRUE)
)

ttl = list(SH = "Agosto", NH = "Enero")
#####  MAPA

lims = list(SH = c(-15,10), NH = c(-5,15))
orien = list(SH = c(-90,0,0), NH = c(90,0,0))
yl = list(SH = c(-90,lats_hs[1]), NH = c(lats_hn[ny],90))
ybrks = list(SH = c(-30,-60), NH = c(60,30))

MAP = list(SH = map_data("world", ylim = yl$SH), NH = map_data("world", ylim = yl$NH))


########## Grafico ############################################################

df = data.frame(
  lons = rep(lons, ny),
  lats = rep(lats$SH, each = nx),
  pv = array(DATA$SH[,,"pv"]*1e6),
  gh = array(DATA$SH[,,"gh"]/9800)
)

GG_SH = ggplot(df) +
  aes(x = lons, y = lats) +
  geom_contour_fill(aes(z = pv), breaks = brks)+
  scale_fill_gradientn(
    colors = pv_cols,
    breaks = brks, 
    labels = labs,
    limits = range(brks)
  ) +
  guides(
    fill = guide_colourstrip(
      title = " PVU", barwidth = 13, barheight = 0.3, title.position = "right", title.hjust = 15,
      label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
    )
  ) +
  geom_contour2(aes(z = gh), breaks = gbrks, size = 0.2) +
  geom_text_contour(aes(z = gh), breaks = 26:32, size = 2.2, stroke = 0.1) +
  geom_map(data = MAP$SH, aes(x = long, y = lat, map_id = region), map = MAP$SH, fill = NA, color = "black", size = 0.1, inherit.aes = FALSE) +
  coord_map("stereographic", orientation = orien$SH, ylim = yl$SH) +
  scale_y_continuous(breaks = ybrks$SH, labels = NULL, name = NULL) +
  scale_x_continuous(breaks = NULL,name = NULL, labels = NULL) +
  geom_vline(xintercept = seq(-120,180,60), size = 0.1, linetype = "dashed", color = "gray5") +
  geom_ribbon(data = df, aes(x = lons, y = lats, ymin = lims$SH[1], ymax = lims$SH[2]), fill = "white") +
  labs(title = ttl$SH) +
  theme_light() +
  theme(
    panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray5"),
    panel.background = element_blank(), axis.ticks=element_blank(),
    legend.position="bottom",plot.title = element_text(hjust = 0, vjust = -10, size = 8),
    legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
    legend.margin = margin(c(0.05,0,0,0), unit = "cm"),
    legend.box.spacing = unit(0, units = "cm"),
    plot.margin = margin(c(-0.5,0,0,0), unit = "cm")
  )


df = data.frame(
  lons = rep(lons, ny),
  lats = rep(lats$NH, each = nx),
  pv = array(DATA$NH[,,"pv"]*1e6),
  gh = array(DATA$NH[,,"gh"]/9800)
)

GG_NH = ggplot(df) +
  aes(x = lons, y = lats) +
  geom_contour_fill(aes(z = pv), breaks = brks)+
  scale_fill_gradientn(
    colors = pv_cols,
    breaks = brks, 
    labels = labs,
    limits = range(brks)
  ) +
  guides(
    fill = guide_colourstrip(
      title = " PVU", barwidth = 13, barheight = 0.3, title.position = "right", title.hjust = 15,
      label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
    )
  ) +
  geom_contour2(aes(z = gh), breaks = gbrks, size = 0.2) +
  geom_text_contour(aes(z = gh), breaks = 26:32, size = 2.2, stroke = 0.1) +
  geom_map(data = MAP$NH, aes(x = long, y = lat, map_id = region), map = MAP$NH, fill = NA, color = "black", size = 0.1, inherit.aes = FALSE) +
  coord_map("stereographic", orientation = orien$NH, ylim = yl$NH) +
  scale_y_continuous(breaks = ybrks$NH, labels = NULL, name = NULL) +
  scale_x_continuous(breaks = NULL,name = NULL, labels = NULL) +
  geom_vline(xintercept = seq(-120,180,60), size = 0.1, linetype = "dashed", color = "gray5") +
  geom_ribbon(data = df, aes(x = lons, y = lats, ymin = lims$NH[1], ymax = lims$NH[2]), fill = "white") +
  labs(title = ttl$NH) +
  theme_light() +
  theme(
    panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray5"),
    panel.background = element_blank(), axis.ticks=element_blank(),
    legend.position="bottom",plot.title = element_text(hjust = 0, vjust = -10, size = 8),
    legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
    legend.margin = margin(c(0.05,0,0,0), unit = "cm"),
    legend.box.spacing = unit(0, units = "cm"),
    plot.margin = margin(c(-0.5,0,0,0), unit = "cm")
  )


figure = ggarrange(plotlist = list(GG_SH, GG_NH), nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom" )

ggsave(
  here("outputs/antecedentes/HSvsHN.png"),
  figure, device = "png", width = 12, height = 6, units = "cm", dpi = 200
)


























