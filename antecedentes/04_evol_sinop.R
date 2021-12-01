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

########## Predefiniendo variables ############################################
vars = c("t", "gh"); nv = length(vars)
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(0, -90, -2.5); ny = length(lats)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)

DATA = array(
  double(nx*ny*nt*nv), 
  dim = c(nx,ny,nt,nv), 
  dimnames = list(lons = lons, lats = lats, days = dates, vars = vars)
)

for (v in 1:nv){
  nc = nc_open(
    here("inputs/era5_plvls", paste0(vars[v], "_climatology.nc"))
  )
  DATA[,,,v] = ncvar_get(
    nc,
    start = c(
      1,
      which(ncvar_get(nc, "lats") == 0),
      which(ncvar_get(nc, "plvl") == 10),
      1
    ), 
    count = c(-1,-1,1,-1)
    )
  
  nc_close(nc); rm(nc)
}


########## Evolucion QUINCENAL ################################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
########## Climatolagias ######################################################
mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
qini = mini + 15


mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
qname = c(
  "1-Ene", "2-Ene", "1-Feb", "2-Feb", "1-Mar", "2-Mar", "1-Abr", "2-Abr", "1-May", "2-May", "1-Jun", "2-Jun",
  "1-Jul", "2-Jul", "1-Ago", "2-Ago", "1-Sep", "2-Sep", "1-Oct", "2-Oct", "1-Nov", "2-Nov", "1-Dec", "2-Dec"
)
nq = length(qname)

imp = seq(1, length.out = 12, by = 2)
par = seq(2, length.out = 12, by = 2)

CLIM = array(double(nx*ny*nq*nv), dim = c(nx,ny,nq,nv), dimnames = list(lons = lons, lats = lats, qname = qname, vars = vars))

for (v in vars) {
  for (i in 1:12){
    CLIM[,,imp[i],v] = apply(DATA[,,mini[i]:(qini[i]-1),v], c(1,2), mean)
  }
  for (p in 1:12){
    CLIM[,,par[p],v] = apply(DATA[,,qini[p]:(mini[p+1]-1),v], c(1,2), mean)
  }
}

########## Acomodando los datos para graficar #################################
d360 = CLIM["180",,,]
d360 = abind(d360, CLIM[which(lons > 180),,,], along = 1)
d360 = abind(d360, CLIM[which(lons <= 180),,,], along = 1)
lons = seq(-180,180,2.5); nx = length(lons)
dimnames(d360) = list(lons = lons, lats = lats, qname = qname, vars = vars)

########## Graficos ###########################################################
tbrks = seq(-85,-15,5)
gbrks = seq(26.5, 32, 0.5)


GG = vector("list", nq); names(GG) = qname

for (q in qname) {
  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    temp = array(d360[,,q,"t"]) - 273,
    geop = array(d360[,,q,"gh"])/9800
  )

GG[[q]] = ggplot(df) +
  aes(x = lons, y = lats) +
  geom_contour_fill(aes(z = temp), breaks = tbrks) +
  scale_fill_gradientn(
    colors = sequential_hcl(length(tbrks)-1, palette = "YlGnBu"),
    limits = range(tbrks),
    breaks = tbrks
  ) +
  guides(
    fill = guide_colourstrip(
      title = "°C", barwidth = 14, barheight = 0.3, title.position = "right", title.hjust = 15,
      label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
    )
  ) +
  geom_contour2(aes(z = geop), breaks = gbrks, size = 0.2) +
  geom_text_contour(aes(z = geop), breaks = 27:32, size = 1.75, stroke = 0.1) +
  geom_map(data = map_data("world"), map = map_data("world"), aes(map_id = region), fill = NA, color = "white", size = 0.1, inherit.aes = FALSE) +
  coord_map("ortho", orientation = c(-90,0,0)) +
  scale_y_continuous(breaks = c(-60,-30), labels = NULL, name = NULL) +
  scale_x_continuous(breaks = seq(-120,180,60), labels = NULL, name = NULL) +
  annotate("text", x = -120, y = -30, label = "30°S", color = "gray60", size = 1.60) +
  annotate("text", x = -120, y = -60, label = "60°S", color = "gray60", size = 1.60) +
  labs(title = q) +
  theme(
    panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.05, linetype = "dashed", color = "gray15"),
    panel.background = element_blank(), axis.ticks=element_blank(),
    legend.position="bottom",plot.title = element_text(hjust = -0.1, vjust = -5, size = 7),
    legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
    plot.margin = margin(c(-0.1,-0.1,-0.1,-0.1), unit = "cm")
  )
}

figure = ggarrange(plotlist = GG, nrow = 4, ncol = 6, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/antecedentes/evol_sinop_quincenal.png"),
  figure, device = "png", width = 16, height = 12, units = "cm", dpi = 200
)


########## Evolucion DECADAL ASON #############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
########## Climatolagias ######################################################
mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
dmini = c(213,244,274,305,335)
dini2 = dmini + 10
dini3 = dmini + 20

dname = c(
  "1-Ago", "2-Ago", "3-Ago", "1-Sep", "2-Sep", "3-Sep", "1-Oct", "2-Oct", "3-Oct", "1-Nov", "2-Nov", "3-Nov"
)
nd = length(dname)

d1 = seq(1, length.out = 4, by = 3)
d2 = seq(2, length.out = 4, by = 3)
d3 = seq(3, length.out = 4, by = 3)

CLIM = array(double(nx*ny*nd*nv), dim = c(nx,ny,nd,nv), dimnames = list(lons = lons, lats = lats, dname = dname, vars = vars))

for (v in vars){
  for (i in 1:4) {
    CLIM[,,d1[i],v] = apply(DATA[,,dmini[i]:(dini2[i]-1),v], c(1,2), mean)
  }
  for (i in 1:4) {
    CLIM[,,d2[i],v] = apply(DATA[,,dini2[i]:(dini3[i]-1),v], c(1,2), mean)
  }
  for (i in 1:4) {
    CLIM[,,d3[i],v] = apply(DATA[,,dini3[i]:(dmini[i+1]-1),v], c(1,2), mean)
  }
}

########## Acomodando los datos para graficar #################################
d360 = CLIM["180",,,]
d360 = abind(d360, CLIM[which(lons > 180),,,], along = 1)
d360 = abind(d360, CLIM[which(lons <= 180),,,], along = 1)
lons = seq(-180,180,2.5); nx = length(lons)
dimnames(d360) = list(lons = lons, lats = lats, dname = dname, vars = vars)

########## Graficos ###########################################################
tbrks = seq(-80,-15,5)
gbrks = seq(26.5, 32, 0.5)


GG = vector("list", nd); names(GG) = dname

for (d in dname) {
  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    temp = array(d360[,,d,"t"]) - 273,
    geop = array(d360[,,d,"gh"])/9800
  )
  
  GG[[d]] = ggplot(df) +
    aes(x = lons, y = lats) +
    geom_contour_fill(aes(z = temp), breaks = tbrks) +
    scale_fill_gradientn(
      colors = sequential_hcl(length(tbrks)-1, palette = "YlGnBu"),
      limits = range(tbrks),
      breaks = tbrks
    ) +
    guides(
      fill = guide_colourstrip(
        title = "°C", barwidth = 14, barheight = 0.3, title.position = "right", title.hjust = 15,
        label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
      )
    ) +
    geom_contour2(aes(z = geop), breaks = gbrks, size = 0.2) +
    geom_text_contour(aes(z = geop), breaks = 27:32, size = 1.75, stroke = 0.1) +
    geom_map(data = map_data("world"), map = map_data("world"), aes(map_id = region), fill = NA, color = "white", size = 0.1, inherit.aes = FALSE) +
    coord_map("ortho", orientation = c(-90,0,0)) +
    scale_y_continuous(breaks = c(-60,-30), labels = NULL, name = NULL) +
    scale_x_continuous(breaks = seq(-120,180,60), labels = NULL, name = NULL) +
    annotate("text", x = -120, y = -30, label = "30°S", color = "gray60", size = 1.60) +
    annotate("text", x = -120, y = -60, label = "60°S", color = "gray60", size = 1.60) +
    labs(title = d) +
    theme(
      panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.05, linetype = "dashed", color = "gray15"),
      panel.background = element_blank(), axis.ticks=element_blank(),
      legend.position="bottom",plot.title = element_text(hjust = -0.1, vjust = -5, size = 7),
      legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
      plot.margin = margin(c(-0.1,-0.1,-0.1,-0.1), unit = "cm")
    )
}

figure = ggarrange(plotlist = GG, nrow = 3, ncol = 4, common.legend = TRUE, legend = "bottom" )

ggsave(
  here("outputs/antecedentes/evol_sinop_decadal.png"),
  figure, device = "png", width = 16, height = 12, units = "cm", dpi = 200
)


########## Evolucion PENTADAL ################################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
########## Climatolagias ######################################################

pini = seq(213, by = 5, length.out = 25)
pname = dates[pini[1:24]]; np = length(pname)
for (i in 1:24){
  mon = substr(pname[i], 1, 2)
  mon = switch(as.numeric(mon)-7, "Ago", "Sep", "Oct", "Nov")
  pname[i] = paste0(substr(pname[i],4,5), "-", mon)
}

CLIM = array(double(nx*ny*np*nv), dim = c(nx,ny,np,nv), dimnames = list(lons = lons, lats = lats, pname = pname, vars = vars))

for (v in vars) {
  for (i in 1:24) {
    CLIM[,,i,v] = apply(DATA[,,pini[i]:(pini[i+1]-1),v], c(1,2), mean)
  }
}

########## Acomodando los datos para graficar #################################
d360 = CLIM["180",,,]
d360 = abind(d360, CLIM[which(lons > 180),,,], along = 1)
d360 = abind(d360, CLIM[which(lons <= 180),,,], along = 1)
lons = seq(-180,180,2.5); nx = length(lons)
dimnames(d360) = list(lons = lons, lats = lats, pname = pname, vars = vars)

########## Graficos ###########################################################
tbrks = seq(-80,-15,5)
gbrks = seq(26.5, 32, 0.5)


GG = vector("list", np); names(GG) = pname

for (p in pname) {
  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    temp = array(d360[,,p,"t"]) - 273,
    geop = array(d360[,,p,"gh"])/9800
  )
  
  GG[[p]] = ggplot(df) +
    aes(x = lons, y = lats) +
    geom_contour_fill(aes(z = temp), breaks = tbrks) +
    scale_fill_gradientn(
      colors = sequential_hcl(length(tbrks)-1, palette = "YlGnBu"),
      limits = range(tbrks),
      breaks = tbrks
    ) +
    guides(
      fill = guide_colourstrip(
        title = "°C", barwidth = 14, barheight = 0.3, title.position = "right", title.hjust = 15,
        label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
      )
    ) +
    geom_contour2(aes(z = geop), breaks = gbrks, size = 0.2) +
    geom_text_contour(aes(z = geop), breaks = 27:32, size = 1.75, stroke = 0.1) +
    geom_map(data = map_data("world"), map = map_data("world"), aes(map_id = region), fill = NA, color = "white", size = 0.1, inherit.aes = FALSE) +
    coord_map("ortho", orientation = c(-90,0,0)) +
    scale_y_continuous(breaks = c(-60,-30), labels = NULL, name = NULL) +
    scale_x_continuous(breaks = seq(-120,180,60), labels = NULL, name = NULL) +
    annotate("text", x = -120, y = -30, label = "30°S", color = "gray60", size = 1.60) +
    annotate("text", x = -120, y = -60, label = "60°S", color = "gray60", size = 1.60) +
    labs(title = p) +
    theme(
      panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.05, linetype = "dashed", color = "gray15"),
      panel.background = element_blank(), axis.ticks=element_blank(),
      legend.position="bottom",plot.title = element_text(hjust = -0.1, vjust = -5, size = 7),
      legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
      plot.margin = margin(c(-0.1,-0.1,-0.1,-0.1), unit = "cm")
    )
}

figure = ggarrange(plotlist = GG, nrow = 4, ncol = 6, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/antecedentes/evol_sinop_pentadal.png"),
  figure, device = "png", width = 16, height = 12, units = "cm", dpi = 200
)

########## Evolucion MENSUAL ################################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
########## Climatolagias ######################################################
mini = c(213,244,274,305,335)
month = c("Ago", "Sep", "Oct", "Nov"); nm = length(month)

CLIM = array(double(nx*ny*nm*nv), dim = c(nx,ny,nm,nv), dimnames = list(lons = lons, lats = lats, month = month, vars = vars))
for (v in vars){
  for (m in 1:nm) {
    CLIM[,,m,v] = apply(DATA[,,mini[m]:(mini[m+1]-1),v], c(1,2), mean)
  }
}
########## Acomodando los datos para graficar #################################
d360 = CLIM["180",,,]
d360 = abind(d360, CLIM[which(lons > 180),,,], along = 1)
d360 = abind(d360, CLIM[which(lons <= 180),,,], along = 1)
lons = seq(-180,180,2.5); nx = length(lons)
dimnames(d360) = list(lons = lons, lats = lats, month = month, vars = vars)

########## Graficos ###########################################################
tbrks = seq(-80,-15,5)
gbrks = seq(28, 31, 0.5)
SH = map_data("world", ylim = c(-90,lats[1]))

GG = vector("list", nm); names(GG) = month


for (m in month) {
  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    temp = array(d360[,,m,"t"]) - 273,
    geop = array(d360[,,m,"gh"])/9800
  )
  
  GG[[m]] = ggplot(df) +
    aes(x = lons, y = lats) +
    geom_contour_fill(aes(z = temp), breaks = tbrks) +
    scale_fill_gradientn(
      colors = sequential_hcl(length(tbrks)-1, palette = "YlGnBu"),
      limits = range(tbrks),
      breaks = tbrks
    ) +
    guides(
      fill = guide_colourstrip(
        title = "°C", barwidth = 14, barheight = 0.3, title.position = "right", title.hjust = 15,
        label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
      )
    ) +
    geom_contour2(aes(z = geop), breaks = gbrks, size = 0.2) +
    geom_text_contour(aes(z = geop), breaks = 28:31, size = 2.2, stroke = 0.1) +
    geom_map(data = SH, aes(x = long, y = lat, map_id = region), map = SH, fill = NA, color = "white", size = 0.1, inherit.aes = FALSE) +
    coord_map("stereographic", orientation = c(-90,0,0), ylim = c(-90,lats[1])) +
    scale_y_continuous(breaks = c(-60,-30), labels = NULL, name = NULL) +
    scale_x_continuous(breaks = NULL,name = NULL, labels = NULL) +
    geom_vline(xintercept = seq(-120,180,60), size = 0.05, linetype = "dashed", color = "white") +
    annotate("text", x = -120, y = -30, label = "30°S", color = "gray90", size = 2) +
    annotate("text", x = -120, y = -60, label = "60°S", color = "gray90", size = 2) +
    geom_ribbon(data = df, aes(x = lons, y = lats, ymin = lats[1], ymax = 25), fill = "white") +
    labs(title = m) +
    theme(
      panel.ontop=TRUE, panel.border = element_blank(), panel.grid = element_line(size = 0.05, linetype = "dashed", color = "gray15"),
      panel.background = element_blank(), axis.ticks=element_blank(),
      legend.position="bottom",plot.title = element_text(hjust = 0.1, vjust = -10, size = 8),
      legend.text = element_text(size = 6), legend.title = element_text(size = 7, face = "bold"),
      legend.margin = margin(c(0.05,0,0,0), unit = "cm"),
      legend.box.spacing = unit(0, units = "cm"),
      plot.margin = margin(c(-0.5,0,-0.0,0), unit = "cm")
    )
}

figure = ggarrange(plotlist = GG, nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom" )

ggsave(
  here("outputs/antecedentes/evol_sinop_mensual.png"),
  figure, device = "png", width = 16, height = 4.5, units = "cm", dpi = 200
)

