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

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(90,-90,-2.5); ny = length(lats)
plvl = c(1000,700,500,300,200,100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
vars = c("t","u"); nv = length(vars)


########## Extraccion de datos ################################################
DATA = array(double(nx*ny*nz*nt*nv), c(nx,ny,nz,nt,nv), dimnames = list(lons = lons, lats = lats, plvl = plvl, dates = dates, vars = vars))

for (v in 1:nv){
  nc = nc_open(
    here("inputs/era5_plvls", paste0(vars[v], "_climatology.nc"))
  )
  DATA[,,,,v] = ncvar_get(nc); nc_close(nc); rm(nc)
}

########## Climatologias ######################################################
mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")


xmean = apply(DATA, c(2,3,4,5), mean)
CLIM_tri = abind(
  apply(xmean[,,c(1:59, 335:365),], c(1,2,4), mean),
  apply(xmean[,,152:243,], c(1,2,4), mean),
  along = 4
)

tri = c("DEF", "JJA"); names(tri) = tri
dimnames(CLIM_tri) = list(lats = lats, plvl = plvl, vars = vars, tri = tri)

########## Graficos ###########################################################
tbrks = seq(-90,30,10); nb = length(tbrks)
tlabs = as.character(tbrks)
#tlabs[which(even(1:nb))] = " "
axis_y_right = list(DEF = element_blank(), JJA = element_text())
axis_y_left = list(DEF = element_text(), JJA = element_blank())
ticks_y_right = list(DEF = element_blank(), JJA = element_line())
ticks_y_left = list(DEF = element_line(), JJA = element_blank())
title_y_right = list(DEF = element_blank(), JJA = element_text(size = 5, color = "gray40", angle = 0, vjust = 1.06, margin = margin(l = -0.1, unit = "cm")))
title_y_left = list(DEF = element_text(size = 5, color = "gray40", angle = 0, vjust = 1.06, margin = margin(r = -0.3, unit = "cm")), JJA = element_blank())
title_hjust = c(0.98,0.02); names(title_hjust) = tri


GG = vector("list", 2); names(GG) = tri

for (t in tri){
df = data.frame(
  lats = rep(lats, nz),
  plvl = rep(plvl, each = ny),
  temp = array(CLIM_tri[,,"t",t]) - 273,
  uwind = array(CLIM_tri[,,"u",t])
)

GG[[t]] = ggplot(df) +
  aes(x = lats, y = plvl) +
  geom_contour_fill(aes(z = temp), breaks = tbrks) +
  geom_contour2(aes(z = temp), breaks = 0, linetype = "solid", size = 0.09, color = "gray30", na.fill = FALSE) +
  scale_fill_gradientn(
    colours = sequential_hcl(nb-1, palette = "YlGnBu", rev = FALSE),
    breaks = tbrks, labels = tlabs, limits = range(tbrks)
  ) +
  guides(
    fill = guide_colorstrip(
      title = "°C", barwidth = 14, barheight = 0.3, title.position = "right", title.hjust = 15,
      label.vjust = 2, ticks = TRUE, ticks.color = "gray50"
      )
  ) +
  geom_contour2(aes(z = uwind), breaks = seq(10,100,10), linetype = "solid", size = 0.1) +
  geom_contour2(aes(z = uwind), breaks = 0, linetype = "solid", size = 0.2) +
  geom_contour2(aes(z = uwind), breaks = seq(-60,-10,10), linetype = "dashed", size = 0.1) +
  geom_text_contour(aes(z = uwind), breaks = seq(-60,100,20), size = 2.00, stroke = 0.1) +
  scale_y_level(
    breaks = c(1000,500,200,100,50,20,10,5,2,1),
    expand = c(0,0), name = "hPa",
    sec.axis =  dup_axis(
      labels = as.character(c(0,5,10,15,20,25,30,35,40,45)),
      name = "km"
    )
  ) +
  scale_x_latitude(breaks = c(60,30,0,-30,-60)
  ) +
  labs(title = tri[t]) +
  geom_hline(yintercept = c(100,10), color = "white", size = 0.1, linetype = "dashed") +
  geom_vline(xintercept = c(60,-60), color = "white", size = 0.1, linetype = "dashed") +
  theme_light()+
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 7, margin = margin(t = 10)),
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 5),
    axis.title = element_text(size = 6, color = "gray10", face = "bold"),
    axis.text = element_text(size = 5),
    axis.text.y.right = axis_y_right[[t]],
    axis.text.y.left = axis_y_left[[t]],
    axis.ticks.y.left = ticks_y_left[[t]],
    axis.ticks.y.right = ticks_y_right[[t]],
    axis.title.y.left = title_y_left[[t]],
    axis.title.y.right = title_y_right[[t]],
    plot.margin = margin(c(-5,0,3,0)),
    legend.margin = margin(c(-10,0,0,0))
  )
}

figure = ggarrange(plotlist = GG, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/antecedentes/plvl_lats.png"),
  figure, device = "png", width = 12, height = 6.75, units = "cm", dpi = 200
)



