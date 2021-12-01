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
require(zoo)            # Tools for time series

source(here("functions","is_leapyear.R"))   # Ask if year is leap
########## Directories ########################################################
##### Guarda en una lista las rutas de los datos
DIR = list(
  ALLDATA = "../../../datos", # Todos los datos de ERA5 en Server Vegeta
  ERA5_I = "../../../datos/ERA5/day", # Datos diarios para los 12 meses del año de T,U,V,Z,PV desde 1000 a 300 hpa entre 1979-2018.
  ERA5_II = "../../../datos/ERA5_aux/nuevos/daily" # Datos diarios para los 12 meses del año de T,U,V,Z,PV desde 200 a 1 hpa entre 1979-2018 y desde 1000 a 1 hpa para 2019
)

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
days = 152:273; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("06","07","08","09") 
years = 1980:2019; na = length(years)

DATA = array(
  double(nx*nt*na), c(nx,nt,na),
  dimnames = list(lons = lons, days = days, years = years)
)



########## Extraccion de datos ################################################
partial = NULL

for (a in years) {
  print(paste0('Extrayendo año ', a))
  for (m in months) {

    nc = nc_open(
      here(DIR$ERA5_II, paste("u", m, a, "day.nc", sep = "_"))
    )
    data_extracted = ncvar_get(
      nc,
      varid = names(nc$var)[2],
      start = c(
        1,
        which(ncvar_get(nc, "latitude") == -60),
        which(ncvar_get(nc, "level") == 1),
        1
      ),
      count = c(-1,1,1,-1)
    )
    nc_close(nc); rm(nc)
    
    partial = abind(partial, data_extracted, along = 2)
  }
  DATA[,,as.character(a)] = partial
  partial = NULL
}
########## Promedio Zonal #####################################################
U_zm = apply(DATA, 2:3, mean)

########## Graficos ###########################################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")

df = data.frame(
    days = rep(days,na),
    uwind = array(U_zm),
    years = rep(years, each = nt)
)


gg = ggplot(df) +
  aes(x = days, y = uwind) +
  geom_line(size = 0.2) +
  scale_y_continuous(
    name = expression(bar(u)~ "[m/s] (60°S, 1 hPa)"),
    expand = c(0,0), breaks = seq(0,100,50), limits = c(-30,130)
  ) +
  scale_x_continuous(
    name = NULL, expand = c(0,0), breaks = mini, labels = mname, limits = range(days)
  )+
  facet_wrap(vars(years), ncol= 8, nrow = 5) +
  theme_light() +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    panel.spacing = unit(0.05, "cm"),
    panel.border = element_rect(size = 0.2),
    panel.grid = element_line(size = 0.3),
    strip.text = element_text(size = 5, color = "black", margin = margin(t = 2, b = 2))
  )

ggsave(
  here("outputs/ftropos","vacillations_1p60S.png"),
  gg, device = "png", width = 16, height = 9, units = "cm", dpi = 200
)



























