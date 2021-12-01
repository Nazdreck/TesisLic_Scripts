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
lats = seq(5,-85,-2.5); ny = length(lats)
plvl = 300
days = 1:92; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-06-01")),6,10)
months = c("06","07","08") 
years = 1979:2019; na = length(years)

DATA = array(
  double(nx*ny*nt*na), c(nx,ny,nt,na),
  dimnames = list(lons = lons, lats = lats, dates = dates, years = years)
)

########## Extraccion de datos ################################################
##### DATOS 1979-2018
### Troposfera
print('Extrayendo datos geopotencial en 300 hPa, periodo: 1979-2018')
nc_tropos = nc_open( # Abre el archivo nc que contiene los datos de la troposfera
  here(DIR$ERA5_I, paste0("era5.geopotential.day.mean.nc"))
)
for (a in as.character(years[1:(na-1)])) {
  print(paste0('Extrayendo año ', a))
  ext = paste0(a, "-06-01")
  
  tropos = ncvar_get( # Extrae los datos de la troposfera
    nc_tropos, # Archivo nc desde donde se van a extraer los datos
    varid = names(nc_tropos$var)[2], # Indica la variable a extraer
    start = c( # Indica donde comenzar a leer
      1,
      which(ncvar_get(nc_tropos, "lat") == lats[1]),
      which(ncvar_get(nc_tropos, "lev") == plvl), # Posicion de la primera longitud, latitud y nivel de presion
      which(as.character(as.Date(ncvar_get(nc_tropos, "time")/24, origin = "1900-01-01")) == ext) # Posicion de la primera fecha
    ),
    count = c( # Indica la cantidad de valores a leer a lo largo de cada dimension
      -1,ny,1, # Lee todos los datos en las dimensiones x,y,z
      nt # Lee un año entero de datos, cambiar por 366 si el año es biciesto
    )
  )
  
  DATA[,,,a] = tropos
}
nc_close(nc_tropos);rm(nc_tropos, tropos)

##### DATOS 2019
partial = NULL

print('Extrayendo datos de geopotencial de 2019')
for (m in months) {
  nc = nc_open(
    here(DIR$ERA5_II, paste("gh",m,"2019", "day.nc", sep = "_"))
  )
  data_extracted = ncvar_get(
    nc,
    varid = names(nc$var)[2],
    start = c(
      1,
      which(ncvar_get(nc, "latitude") == lats[1]),
      which(ncvar_get(nc, "level") == plvl),
      1
    ),
    count = c(
      -1,ny,1,-1
    )
  )
  nc_close(nc); rm(nc)
  
  partial = abind(partial, data_extracted, along = 3)
}

DATA[,,,"2019"] = partial

rm(data_extracted, partial)

########## Analisis Espectral #################################################
FOURIER = array(
  double(ny*nt*na), c(ny,nt,na),
  dimnames = list(lats = lats, dates = dates, years = years)
)
for (y in 1:ny) {
  for (t in 1:nt) {
    for (a in 1:na) {
      FOURIER[y,t,a] = FitWave(DATA[,y,t,a]/9.8, k = 1)$amplitude
    }
  }
}

FOURIER = apply(FOURIER, c(1,3), mean)

########## Percentiles ########################################################
pers = c(10,30,70,90); np = length(pers)

PER = array(double(np*ny), c(ny,np), dimnames = list(lats = lats,percentile = pers))

for (d in 1:ny) {
  for (p in pers) {
    x1 = which(Percentile(FOURIER[d,]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
    x2 = min(FOURIER[d,][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
    PER[d,as.character(p)] = x2 # Guarda el valor minimo encontrado
  }
}

########## Graficos ###########################################################
NINO = as.character(c(1982,1987,1991,1997,2002,2004,2009,2015))
NINA = as.character(c(1985,1988,1998,1999,2000,2007,2010,2011))

df_lines = data.frame(
  lats = rep(lats,11),
  vars = c(
    p10 = PER[,"10"],
    p90 = PER[,"90"],
    p30 = PER[,"30"],
    p70 = PER[,"70"],
    max = apply(FOURIER, 1, max),
    min = apply(FOURIER, 1, min),
    clim = apply(FOURIER, 1, mean),
    d2002 = FOURIER[,"2002"],
    d2019 = FOURIER[,"2019"],
    dnino = apply(FOURIER[,NINO], 1, mean),
    dnina = apply(FOURIER[,NINA], 1, mean)
  ),
  group = rep(LETTERS[1:11], each = ny),
  cols = c(
    rep("10-90%", each = ny*2), rep("30-70%", each = ny*2), rep("Máx/mín", each = ny*2),
    rep("Clim.", each = ny), rep("2002", each = ny), rep("2019", each = ny),
    rep("El Niño", each = ny), rep("La Niña", each = ny)
  )
)

df_lines$cols = factor(df_lines$cols, levels = c("10-90%", "30-70%", "Máx/mín", "Clim.", "El Niño", "La Niña","2002", "2019"))

df_ribbons = data.frame(
  lats = rep(lats, 2),
  pmin = c(p10 = PER[,"10"], p30 = PER[,"30"]),
  pmax = c(p90 = PER[,"90"], p70 = PER[,"70"]),
  group = rep(LETTERS[1:2], each = ny)
)

GG = ggplot() +
  geom_ribbon(
    data = df_ribbons, aes(x = lats, ymin = pmin, ymax = pmax, fill = group), 
    alpha = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("A"= "gray80", "B" = "gray40")) +
  geom_line(
    data = df_lines, aes(x = lats, y = vars, group = group, color = cols, size = cols)
  ) +
  scale_color_manual(
    name = "Leyenda",
    values = c(
      "Máx/mín" = "darkgreen", "Clim." = "black", 
      "2019" = "blue", "2002"= "red", 
      "10-90%"= "gray70", "30-70%" = "gray30",
      "El Niño" = "orange", "La Niña" = "purple"
    )
  ) +
  scale_size_manual(
    values = c(
      "Máx/mín" = 0.15, "Clim." = 0.3, 
      "2019" = 0.3, "2002"= 0.3, 
      "10-90%"= 0.1, "30-70%" = 0.1,
      "El Niño" = 0.2, "La Niña" = 0.2
    )
  ) +
  guides(size = "none") +
  scale_y_continuous(name = NULL, breaks = seq(0,200,50), limits = c(0,210), expand = c(0,0))+
  scale_x_latitude(name = NULL, breaks = seq(-75,-15,15), limits = c(-85,-5)) +
  labs(
    title = expression("Amplitud Z"[1]~"[m] (300 hPa, JJA)")
  ) +
  theme_light() +
  theme(
    plot.margin = margin(0,0,0,0,"cm"),
    plot.title = element_text(size = 6, vjust = -7, hjust = 0.02, face = "bold", margin = margin(0,0,0,0, "cm")),
    legend.key.size = unit(0.275, "cm"),
    legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
    legend.margin = margin(-2, 3, 3, 3),
    legend.key = element_rect(fill = alpha("orange", 0.05)),
    legend.title = element_blank(),
    legend.position = c(0.99,0.97),
    legend.justification = c("right", "top"),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
  )

ggsave(
  here("outputs/ftropos", "Z1_series.png"),
  GG, device = "png", width = 12, height = 5.0625, units = "cm", dpi = 200
)




















