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

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(-35,-85,-2.5); ny = length(lats)
plvl = 10  
days = 91:304; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = c(2002,2019); na = length(years)
vars = c("t","u","v","pv"); nv = length(vars)

DATA = array(
  double(nx*ny*nt*na*nv), c(nx,ny,nt,na,nv),
  dimnames = list(lons = lons, lats = lats, dates = dates, years = years, vars = vars)
)

########## Extrayendo datos ###################################################
for (v in vars) {
  for (a in as.character(years)) {
    nc = nc_open(
      here("inputs/era5_plvls", paste0(v, "_", a, ".nc"))
    )
    
    DATA[,,,a,v] = ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == lats[1]),
        which(ncvar_get(nc, "plvl") == plvl[1]),
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == dates[1])
      ),
      count = c(
        c(-1,ny,1,nt)
      )
    )
    nc_close(nc); rm(nc)
  }
}

# Carga los dias con minimos de flujo de calor
load("inputs/heatflux_mindays.RData")
########## Flujos meridional de VP ############################################
QV = array(
  double(ny*nt*na), c(ny,nt,na),
  dimnames = list(lats = lats, dates = dates,years = years)
)

for (a in 1:na) {
  for (t in 1:nt) {
    for (y in 1:ny) {
      QV[y,t,a] = cov(DATA[,y,t,a,"v"],DATA[,y,t,a,"pv"])
    }
  }
}

##### calculo largo
Q = apply(DATA[,,,,"pv"], c(2,3,4), mean)
V = apply(DATA[,,,,"v"], c(2,3,4), mean)

qprime = DATA[,,,,"pv"]*0
vprime = DATA[,,,,"v"]*0
for (x in 1:nx) {
qprime[x,,,] = DATA[x,,,,"pv"] - Q
vprime[x,,,] = DATA[x,,,,"v"] - V
}

qv = qprime*vprime
QV = apply(qv, c(2,3,4), mean)
########## Medias zonales de T y U ############################################
TK = apply(DATA[,,,,"t"], c(2,3,4), mean)
U = apply(DATA[,,,,"u"], c(2,3,4), mean)

########## Elementos para graficar ############################################
mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")

brks = seq(-3,5.5,0.5); nb = length(brks)
labs = as.character(brks)
labs[which(even(1:nb))] = " "

brks_t = seq(190,260,15)

nn = length(which(brks < 0))
np = length(which(brks > 0))

cols = c(
  sequential_hcl(nn, palette = "YlGn", rev = FALSE),
  sequential_hcl(np, palette = "YlOrRd", rev = TRUE)
)

hfmin = list("2002" = mins2002, "2019" = mins2019)
##### theming
# x axis
tm.axis_text_x_bottom = list("2002" = element_text(), "2019" = element_blank())
tm.axis_ticks_x_bottom = list("2002" = element_line(), "2019" = element_blank())
tm.axis_ticks_x_top = list("2002" = element_blank(), "2019" = element_line())

########## Grafico ############################################################
GG = vector("list", na); names(GG) = years

for (a in as.character(years)) {
  df = data.frame(
    lats = rep(lats, nt),
    days = rep(days, each = ny),
    qv = array(QV[,,a]*1000),
    uwind = array(U[,,a]),
    temp = array(TK[,,a])
  )
  
  GG[[a]] = ggplot(df) +
    aes(x = days, y = lats) +
    geom_contour_fill(aes(z = qv), breaks = brks) +
    scale_fill_gradientn(
      colors = cols, 
      breaks = brks,
      labels = labs,
      limits = range(brks)
    ) +
    guides(
      fill = guide_colorstrip(
        title = expression(bar("q'v'")~ "[10"^3~"K.m"^3~".kg"^{-1}~".s"^{-2}~"]    "), 
        ticks = TRUE, ticks.color = "gray30", barwidth = 10, barheight = 0.3, label.vjust = 2
      )
    ) +
    geom_contour2(aes(z = uwind), breaks = 0, linetype = "solid", size = .3, alpha = 0.6, color = "black") +
    geom_contour2(aes(z = uwind), breaks = seq(20,100,20), linetype = "solid", size = .15, alpha = 0.6, color = "gray20") +
    geom_contour2(aes(z = uwind), breaks = seq(-40,-20,20), linetype = "dashed", size = .15, alpha = 0.6, color = "gray20") +
    geom_contour2(aes(z = temp), breaks = brks_t, linetype = "solid", size = 0.25, alpha = 0.6, color = "mediumblue")+
    geom_text_contour(aes(z = temp), breaks = brks_t, size = 2.2, stroke = 0.1) +
    scale_x_continuous(
      breaks = mini, labels = mname, limits = c(106,289),
      expand = c(0,0), name = NULL, sec.axis = dup_axis(name = NULL)
    ) +
    scale_y_latitude(
      breaks = seq(-75,-45,10), limits = c(-85,-35)
    ) +
    geom_vline(xintercept = hfmin[[a]], color = "darkgreen", linetype = "dashed", size = 0.15) +
    annotate("text", x = 113, y = -40, label = a, color = "black", size = 2) +
    theme_light() +
    theme(
      panel.ontop = TRUE, 
      panel.background = element_blank(), 
      panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray50"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0.05, r = 0.1, b = 0, l = 0, unit = "cm"),
      plot.title =  element_blank(),
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
      # axis.text.y.left = tm.axis_text_y_left[[a]],
      # axis.text.y.right = tm.axis_text_y_right[[a]],
      # axis.ticks.y.left = tm.axis_ticks_y_left[[a]],
      # axis.ticks.y.right = tm.axis_ticks_y_right[[a]],
      axis.text = element_text(size = 5)
    )
  
}

figure = ggarrange(plotlist = GG, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/ftropos/QV_plvl-lats_field.png"),
  figure, device = "png", width = 16, height = 9, units = "cm", dpi = 200
)




















