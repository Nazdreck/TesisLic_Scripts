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
plvl =c(1000,700,500,300,200,100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
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
CLIM_10 = array(double(nt*ny*nv), c(nt,ny,nv), dimnames = list(dates = dates, lats = lats, vars = vars))
CLIM_10[,"0",] = NA
CLIM_10[1:184,which(lats > 0),] = aperm(xmean[which(lats > 0), "10", 182:365,], c(2,1,3))
CLIM_10[185:365,which(lats > 0),] = aperm(xmean[which(lats > 0), "10", 1:181,], c(2,1,3))
CLIM_10[,which(lats < 0),] = aperm(xmean[which(lats < 0), "10",,], c(2,1,3))


CLIM_60 = abind(
  aperm(xmean["-60",,,], c(2,1,3)),
  abind(
    aperm(xmean["60",,182:365,], c(2,1,3)),
    aperm(xmean["60",,1:181,], c(2,1,3)),
    along = 1
  ),
  along = 4
)

hem = c("S", "N")
dimnames(CLIM_60) = list(dates = dates, plvl = plvl, vars = vars, hem = hem)

########## Graficos ###########################################################

########## Grafico latitud vs tiempo
ubrks = seq(-90,90,10); nb = length(ubrks)
ulabs = as.character(ubrks)
ulabs[which(even(1:nb))] = " "

colors = c(
  sequential_hcl((nb-1)/2, palette = "PuRd", rev = FALSE), # YlGnBu / Blues /
  sequential_hcl((nb-1)/2, palette = "Greens", rev = TRUE) # YlOrRd / Reds / OrRd
)

df = data.frame(
  days = rep(days, ny),
  lats = rep(lats, each = nt),
  temp = array(CLIM_10[,,"t"])-273,
  uwind = array(CLIM_10[,,"u"])
)

GG_lats = ggplot(df) +
  aes(x = days, y = lats) +
  geom_contour_fill(aes(z = uwind), breaks = ubrks, na.fill = FALSE) +
  geom_contour2(aes(z = uwind), breaks = 0, linetype = "dotdash", size = 0.15, color = "darkgreen", na.fill = FALSE) +
  scale_fill_gradientn(
    #colours = sequential_hcl(nb-1, palette = "viridis", rev = FALSE),
    colors = colors,
    breaks = ubrks, labels = ulabs, limits = range(ubrks)
  ) +
  guides(
    fill = guide_colorstrip(
      title = "m/s   ", barwidth = 14, barheight = 0.3, title.position = "left", title.hjust = 0.5,
      ticks = TRUE, ticks.color = "gray50", label.vjust = 2
    )
  ) +
  geom_contour2(aes(z = temp), breaks = seq(-90,-10,10), linetype = "dashed", size = 0.1, na.fill = FALSE) +  
  geom_text_contour(aes(z = temp), breaks = seq(-90,-10,10), size = 2.2, stroke = 0.1) +
  geom_ribbon(aes(ymin = -5, ymax = 5), fill = "gray75") +
  scale_x_continuous(
    breaks = mini, labels = mname,
    expand = c(0,0), name = NULL,
    sec.axis =  dup_axis(
      labels = c(mname[7:12], mname[1:7]),
      name = NULL
    )
  ) +
  scale_y_latitude(
    breaks = c(90,60,30,0,-30,-60,-90)
  ) +
  geom_hline(
    yintercept = c(-60,60), color = "black", linetype = "dashed", size = 0.1
  ) +
  geom_vline(
    xintercept = c(80,264), color = "black", linetype = "dashed", size = 0.1
  ) +
  theme_light() +
  theme(
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    axis.text = element_text(size = 5),
    plot.margin = margin(l = -9, r = 5)
  )

########## Grafico presion vs tiempo
ax_ticks_xtop = list(S = element_line(), N = element_blank())
ax_ticks_xbot = list(S = element_blank(), N = element_line())

GG_plvl = vector("list", 2); names(GG_plvl) = hem

for (h in hem) {
  df = data.frame(
    days = rep(days, nz),
    plvl = rep(plvl, each = nt),
    temp = array(CLIM_60[,,"t", h])-273,
    uwind = array(CLIM_60[,,"u",h])
  )
  
  GG_plvl[[h]] = ggplot(df) +
    aes(x = days, y = plvl) +
    geom_contour_fill(aes(z = uwind), breaks = ubrks) +
    geom_contour2(aes(z = uwind), breaks = 0, linetype = "dotdash", size = 0.15, color = "darkgreen", na.fill = FALSE) +
    scale_fill_gradientn(
      #colours = sequential_hcl(nb-1, palette = "viridis", rev = FALSE),
      colors = colors,
      breaks = ubrks, labels = ulabs, limits = range(ubrks)
    ) +
    guides(
      fill = guide_colorstrip(
        title = "m/s   ", barwidth = 14, barheight = 0.3, title.position = "left", title.hjust = 0.5,
        ticks = TRUE, ticks.color = "gray50", label.vjust = 2
      )
    ) +
    # geom_contour2(aes(z = temp), breaks = seq(10,20,10), linetype = "solid", size = 0.1, na.fill = FALSE) +
    geom_contour2(aes(z = temp), breaks = 0, linetype = "solid", size = 0.2, na.fill = FALSE) +
    geom_contour2(aes(z = temp), breaks = seq(-90,-10,10), linetype = "dashed", size = 0.1, na.fill = FALSE) +  
    geom_text_contour(aes(z = temp), breaks = seq(-90,0,10), size = 2.2, stroke = 0.1) +
    scale_y_level(
      breaks = c(1000,500,200,100,50,20,10,5,2,1), name = "hPa"
    ) +
    scale_x_continuous(name = NULL, limits = range(1:365), breaks = mini, labels = mname, expand = c(0,0), sec.axis = dup_axis(name = NULL)) +
    geom_hline(
      yintercept = c(100,10), color = "black", linetype = "dashed", size = 0.1
    ) +
    geom_vline(
      xintercept = c(80,264), color = "black", linetype = "dashed", size = 0.1
    ) +
    theme_light() +
    theme(
      legend.title = element_text(face = "bold", size = 7),
      legend.text = element_text(size = 5),
      legend.position = "bottom",
      axis.title.y = element_text(size = 5, color = "gray40", angle = 0, margin = margin(r = -0.6, unit = "cm")),
      axis.text = element_text(size = 5),
      axis.text.x = element_blank(),
      axis.ticks.x.top = ax_ticks_xtop[[h]],
      axis.ticks.x.bottom = ax_ticks_xbot[[h]],
      plot.margin = margin(l = -9, r = 5)
    )
}


figure = ggarrange(plotlist = list(GG_plvl[[2]], GG_lats, GG_plvl[[1]]), ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom",align = "v")

ggsave(
  here("outputs/antecedentes/evol_temp.png"),
  figure, device = "png", width = 9, height = 16, units = "cm", dpi = 200
)


