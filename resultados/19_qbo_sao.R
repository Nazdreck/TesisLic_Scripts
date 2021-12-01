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
require(gtools)

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(5,-5,-2.5); ny = length(lats)
plvl = c(200,100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
days = 0:(365*20+4); nt = length(days)
dates = as.character(as.Date(days, origin = "2000-01-01"))
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = 2000:2019; na = length(years)

########## Extraccion de datos ################################################
#DATA = array(double(nx*ny*nz*nt), c(nx,ny,nz,nt), dimnames = list(lons = lons, lats = lats, plvl = plvl, dates = dates))
DATA = array(double(nz*nt), c(nt, nz), dimnames = list(dates = dates, plvl = plvl))

DIR = "../../../datos/ERA5_aux/nuevos/daily"

for (a in years) {
  for (m in months) {
    nc = nc_open(
      here(DIR, paste("u", m, a, "day.nc", sep = "_"))
    )
    data_extracted = ncvar_get(
      nc,
      varid = names(nc$var)[2],
      start = c(
        1,
        which(ncvar_get(nc, "latitude") == 5),
        which(ncvar_get(nc, "level") == 200),
        1
      ),
      count = c(-1,5,-1,-1)
    )
    date_extracted = as.character(as.Date(ncvar_get(nc, "time")/24, origin = "1900-01-01"))
    nc_close(nc); rm(nc)
    
    DATA[date_extracted,] = t(apply(data_extracted, c(3,4), mean))
  }
}

########## Media cada x dias ##################################################
x = 5
pts = seq(from = 1, by = x, length.out = nt/x); np = length(pts)
MEAN = array(double(nz*np), c(np,nz), dimnames = list(time = 1:np, plvl = plvl))

for (i in 1:np) {
  MEAN[which(pts == pts[i]),] = apply(DATA[pts[i]:(pts[i]+4),], c(2), mean)
}
  
########## Grafico ############################################################
brks = seq(-90,90,5); nb = length(brks)
labs = as.character(brks)
labs[which(even(1:nb))] = " "

df = data.frame(
  time = rep(1:np, nz),
  plvl = rep(plvl, each = np),
  uwind = array(MEAN, np*nz)
)
  
gg = ggplot(df) +
  aes(x = time, y = plvl, z = uwind) +
  geom_contour_fill(
    breaks = brks
  ) +
  scale_fill_gradientn(
    colours = divergingx_hcl(nb-1, palette = "PuOr", rev = TRUE),
    breaks = brks, labels = labs, limits = range(brks)
  ) +
  guides(
    fill = guide_colorstrip(
      title = "m/s   ", ticks = TRUE, ticks.color = "gray50", barwidth = 14, barheight = 0.3,
      title.position = "left", label.vjust = 2, title.vjust = 1.15
    )
  ) +
  geom_contour2(breaks = 0, linetype = "solid", size = 0.2) +
  geom_contour2(breaks = seq(10,90,30), linetype = "solid", size = 0.1) +
  geom_contour2(breaks = seq(-90,-10,30), linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = 1000/5, color = "red", size = 0.3, linetype = "solid") +
  geom_vline(xintercept = 7200/5, color = "red", size = 0.3, linetype = "solid") +
  scale_y_level(
    breaks = c(200,100,70,50,30,20,10,7,5,3,2,1), limits = c(100,1), name = "hPa"
  ) +
  scale_x_continuous(
    breaks = seq(1, by = 73, length.out = 21),
    labels = 2000:2020,
    limits = c(1,1461),
    expand = c(0,0), name = NULL
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 5),
    axis.title.y.left = element_text(size = 6, color = "gray50", face = "bold", angle = 0, vjust = 0.98, hjust = 1, margin = margin(r = -0.15, unit = "cm")),
    axis.text = element_text(size = 5),
    plot.margin = margin(c(2,6,0,0)),
    legend.margin = margin(c(-6,0,0,0))
  )

ggsave(
  here("outputs/ftropos/qbo_sao2.png"),
  gg, device = "png", width = 16, height = 6.75, units = "cm", dpi = 200
)
  
