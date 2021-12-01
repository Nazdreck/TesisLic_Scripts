########## Limpia el espacio de trabajo #######################################
rm(list=ls())

########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ncdf4)          # For manipulate .nc files
require(ggplot2)        # Data Visualisation
require(abind)          # Combine Multidimensional Arrays
require(ggpubr)         # "ggplot2" Based Publication Ready Plots
require(metR)           # Tools for Meteorological Analysis
require(scales)         # Scale Functions for visualization
require(colorspace)     # Paletas de Colores.
require(zoo)            # Tools for time series
require(gtools)

source(here("functions","is_leapyear.R"))   # Ask if year is leap
source(here("functions","QGEPflux.R"))   # Ask if year is leap

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(-5,-85,-2.5); ny = length(lats)
plvl = c(100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
days = 1:122; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-06-01")),6,10)
months = c("06","07","08","09") 
years = c(2002,2019); na = length(years)
vars = c("t","u","v"); nv = length(vars)

DATA = array(
  double(nx*ny*nz*nt*na*nv), c(nx,ny,nz,nt,na,nv),
  dimnames = list(lons = lons, lats = lats, plvl = plvl, days = days, years = years, vars = vars)
)

########## Extrayendo datos ###################################################
for (v in vars) {
  for (a in as.character(years)) {
    nc = nc_open(
      here("inputs/era5_plvls", paste0(v, "_", a, ".nc"))
    )
    
    DATA[,,,,a,v] = ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == lats[1]),
        which(ncvar_get(nc, "plvl") == plvl[1]),
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == dates[1])
      ),
      count = c(
        c(-1,ny,-1,nt)
      )
    )
    nc_close(nc); rm(nc)
  }
}

########## Flujos de Calor ####################################################
HF = array(
  double(ny*nz*nt*na), c(ny,nz,nt,na),
  dimnames = list(lats = lats, plvl = plvl, days = days,years = years)
)

for (a in 1:na) {
  for (t in 1:nt) {
    for (z in 1:nz) {
      for (y in 1:ny) {
        HF[y,z,t,a] = cov(DATA[,y,z,t,a,"t"],DATA[,y,z,t,a,"v"])
      }
    }
  }
}
########## Divergencia del Flujo de Eliassen-Palm #############################
divF = array(
  double(ny*nz*nt*na), c(ny,nz,nt,na), 
  dimnames = list(lats = lats, plvl = plvl, dates = dates, years = years)
)

for (a in 1:na){ 
  for (t in 1:nt) {
    divF[,,t,a] = QGEPflux(dataset = DATA[,,,t,a,], longitudes = lons, latitudes = lats, pressure = plvl)[,,4]
  }
}

########## Promedios ##########################################################
proms = c("JJA", "Sep"); np = length(proms)
dds = list(1:92, 93:nt)

HFP = array(
  double(ny*nz*np*na), c(ny,nz,np,na),
  dimnames = list(lats = lats, plvl = plvl, proms = proms,years = years)
)
divFP = array(
  double(ny*nz*np*na), c(ny,nz,np,na),
  dimnames = list(lats = lats, plvl = plvl, proms = proms,years = years)
)
UP = array(
  double(ny*nz*np*na), c(ny,nz,np,na),
  dimnames = list(lats = lats, plvl = plvl, proms = proms,years = years)
)

HFP[,,"JJA",] = apply(HF[,,dds[[1]],], c(1,2,4), mean)
HFP[,,"Sep",] = apply(HF[,,dds[[2]],], c(1,2,4), mean)
divFP[,,"JJA",] = apply(divF[,,dds[[1]],], c(1,2,4), mean)
divFP[,,"Sep",] = apply(divF[,,dds[[2]],], c(1,2,4), mean)
UP[,,"JJA",] = apply(DATA[,,,dds[[1]],,"u"], c(2,3,5), mean)
UP[,,"Sep",] = apply(DATA[,,,dds[[2]],,"u"], c(2,3,5), mean)

########## Elementos para graficar ############################################
brks_vt = seq(-220,-20,20); nb = length(brks_vt)
labs = as.character(brks_vt)
labs[which(odd(1:nb))] = " "

brks_df = seq(-60,-10,10)

cols = sequential_hcl(nb-1, palette = "YlOrRd")

ttl = c("JJA" = "Jun-Jul-Ago", "Sep" = "Septiembre")

A = matrix(
  c(rep(c("2002","2019"),each=2), rep(proms,2)), 
  ncol = 2, nrow = 4
)

##### theming
tm.plot_title = list("2002" = element_text(size = 6, face = "bold", vjust = -8, hjust = 0.03, margin = margin(c(-0.19,0,0,0), unit = "cm")), "2019" = element_blank())
# x axis
tm.axis_text_x_bottom = list("2002" = element_text(), "2019" = element_blank())
tm.axis_ticks_x_bottom = list("2002" = element_line(), "2019" = element_blank())
tm.axis_ticks_x_top = list("2002" = element_blank(), "2019" = element_line())
# y axis
tm.axis_text_y_left = list("JJA"  = element_text(), "Sep" = element_blank())
tm.axis_text_y_right = list("JJA"  = element_blank(), "Sep" = element_text())
tm.axis_ticks_y_left = list("JJA"  = element_line(), "Sep" = element_blank())
tm.axis_ticks_y_right = list("JJA"  = element_blank(), "Sep" = element_line())
########## Grafico ############################################################
GG = vector("list", np*na)
for (i in 1:(np*na)) {
  a = as.character(A[i,1])
  p = A[i,2]
  
  df = data.frame(
    lats = rep(lats, nz),
    plvl = rep(plvl, each = ny),
    hf = array(HFP[,,p,a]),
    dfep = array(divFP[,,p,a]*86400),
    uwind = array(UP[,,p,a])
  )
  
  GG[[i]] = ggplot(df) +
    aes(x = lats, y = plvl) +
    geom_contour_fill(aes(z = hf), breaks = brks_vt) +
    scale_fill_gradientn(
      colors = cols, na.value = "white",
      breaks = brks_vt,
      labels = labs,
      limits = range(brks_vt)
    ) +
    guides(
      fill = guide_colorstrip(
        title = expression(bar("v'T'")~ "[K ms"^{-1}*"]   "), ticks = TRUE, ticks.color = "gray30", barwidth = 10, barheight = 0.3, label.vjust = 2
      )
    ) +
    geom_contour2(aes(z = uwind), breaks = 0, linetype = "solid", size = .3, alpha = 0.6, color = "black") +
    geom_contour2(aes(z = uwind), breaks = seq(10,100,10), linetype = "solid", size = .15, alpha = 0.6, color = "gray20") +
    geom_contour2(aes(z = uwind), breaks = seq(-40,-10,10), linetype = "dashed", size = .15, alpha = 0.6, color = "gray20") +
    geom_contour2(aes(z = dfep), breaks = brks_df, linetype = "dotdash", size = 0.2, color = "mediumblue")+
    geom_text_contour(aes(z = dfep), breaks = brks_df, size = 2.2, stroke = 0.1) +
    scale_y_level(
      name = NULL, breaks = c(100,70,50,30,20,10,7,5,3,2,1), sec.axis = dup_axis(name = NULL)
    ) +
    scale_x_latitude(
      name = NULL, breaks = seq(-80,-20,20), sec.axis = dup_axis(name = NULL)
    ) +
    annotate("text", x = -80, y = 81, label = a, color = "black", size = 2) +
    labs(title = ttl[p]) +
    theme_light() +
    theme(
      panel.ontop = TRUE, 
      panel.background = element_blank(), 
      panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray50"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      plot.title =  tm.plot_title[[a]],
      legend.title = element_text(size = 5, face = "bold"),
      legend.text = element_text(size = 4),
      legend.position = "bottom",
      legend.margin = margin(c(0,0,0,0), unit = "cm"),
      legend.box.spacing = unit(0, units = "cm"),
      axis.title = element_blank(),
      # x-axis
      axis.text.x.bottom = tm.axis_text_x_bottom[[a]],
      axis.text.x.top = element_blank(),
      axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[a]],
      axis.ticks.x.top = tm.axis_ticks_x_top[[a]],
      # y-axis
      axis.text.y.left = tm.axis_text_y_left[[p]],
      axis.text.y.right = tm.axis_text_y_right[[p]],
      axis.ticks.y.left = tm.axis_ticks_y_left[[p]],
      axis.ticks.y.right = tm.axis_ticks_y_right[[p]],
      axis.text = element_text(size = 5)
    )
} 


figure = ggarrange(plotlist = GG, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/ftropos/HF_plvl-lats_field3.png"),
  figure, device = "png", width = 12, height = 9, units = "cm", dpi = 200
)
