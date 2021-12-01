########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ncdf4)          # For manipulate .nc files
require(ggplot2)        # Data Visualisation
require(metR)           # Tools for Meteorological Analysis
require(colorspace)     # Paletas de Colores.
require(abind)          # Combine Multidimensional Arrays
require(ggpubr)         # "ggplot2" Based Publication Ready Plots
require(gtools)

source(here("functions","is_leapyear.R"))   # Ask if year is leap
source(here("functions","QGEPflux.R"))   # Ask if year is leap


########## Directories ########################################################
##### Guarda en una lista las rutas de los datos
DIR = list(
  ALLDATA = "../../../datos", # Todos los datos de ERA5 en Server Vegeta
  ERA5_I = "../../../datos/ERA5/day", # Datos diarios para los 12 meses del año de T,U,V,Z,PV desde 1000 a 300 hpa entre 1979-2018.
  ERA5_II = "../../../datos/ERA5_aux/nuevos/daily" # Datos diarios para los 12 meses del año de T,U,V,Z,PV desde 200 a 1 hpa entre 1979-2018 y desde 1000 a 1 hpa para 2019
)

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(-45,-75,-2.5); ny = length(lats)
plvl = c(1000,850,700,500,400,300,200,100); nz = length(plvl)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("01","02","03","04","05","06","07","08","09","10","11","12") 
vars = c("t","v"); nv = length(vars)
years = c(2002,2019); na = length(years)

DATA = array(
  double(nx*ny*nz*nt*nv*na), c(nx,ny,nz,nt,nv,na),
  dimnames = list(lons = lons, lats = lats, plvl = plvl, dates = dates, vars = vars, years = years)
)

########## Archivos de Datos ##################################################
##### Lista de nombres de variables 
# Primera columna es para ERA5_I. Segunda columna es para era "ERA5_II"
variables = list( 
  Temp = c("temperature","t"), 
  Vwind = c("v_component_of_wind", "v")
)

########## Extraccion de datos ################################################
##### DATOS 2002
### Extraccion de datos de la troposfera
start_date = "2002-01-01" # Escribe en formato fecha el primer dia del año a extraer

for (v in 1:nv) {
  print(paste0("Extrayendo datos para ", vars[v]))
##### DATOS 2002
  print('Estrayendo datos de 2002')
  print("Extrayendo datos de la troposfera")
  
  #print("Abriendo archivo .nc")
  nc_tropos = nc_open( # Abre el archivo nc que contiene los datos de la troposfera
    here(DIR$ERA5_I, paste0("era5.", variables[[v]][1], ".day.mean.nc"))
  )
  
  #print(paste0("Extrayendo 1 año de datos desde el ", start_date, " ..."))
  tropos = ncvar_get( # Extrae los datos de la troposfera
    nc_tropos, # Archivo nc desde donde se van a extraer los datos
    varid = names(nc_tropos$var)[2], # Indica la variable a extraer
    start = c( # Indica donde comenzar a leer
      1,
      which(ncvar_get(nc_tropos, "lat") == lats[1]),
      1, # Posicion de la primera longitud, latitud y nivel de presion
      which(as.character(as.Date(ncvar_get(nc_tropos, "time")/24, origin = "1900-01-01")) == start_date) # Posicion de la primera fecha
    ),
    count = c( # Indica la cantidad de valores a leer a lo largo de cada dimension
      -1,ny,-1, # Lee todos los datos en las dimensiones x,y,z
      365 # Lee un año entero de datos, cambiar por 366 si el año es biciesto
    )
  )
  
  print("Datos de troposfera extraidos")
  
  nc_close(nc_tropos) # Cierra el archivo nc con los datos de troposfera
  rm(nc_tropos) # Remueve la variable nc_tropos del espacio de trabajo
  
  ### Extraccion de datos de la estratosfera
  stratos = NULL # variable donde se van a guardar los datos de estratosfera
  
  print("Extrayendo datos de la estratosfera")
  
  for (m in months) { # Para cada mes
    #print(paste0("Abriendo archivo .nc para el mes ", m))
    nc_stratos = nc_open( # Abre el archivo nc que contiene los datos de la estratosfera
      here(DIR$ERA5_II, paste(variables[[v]][2], m, "2002_day.nc", sep = "_"))
    )
    
    #print("Extrayendo datos ...")
    data_extracted = ncvar_get( # Extrae los datos de la estratosfera
      nc_stratos, # Archivo nc desde donde se van a extraer los datos
      varid = names(nc_stratos$var)[2], # Indica la variable a extraer
      start = c(
        1,
        which(ncvar_get(nc_stratos, "latitude") == lats[1]),
        1,
        1
        ), # Posicion de la primera longitud, latitud y nivel de presion
      count = c(-1,ny,2,-1) # Indica la cantidad de valores a leer a lo largo de cada dimension. Extrae solo 200 y 100 hPa
    )
    nc_close(nc_stratos) # Cierra el archivo nc con los datos de estratosfera
    rm(nc_stratos) # Remueve la variable nc_stratos del espacio de trabajo
    
    #print(paste0("Datos del mes ", m, " extraidos. Combinando..."))
    stratos = abind(stratos, data_extracted, along = 4) # Pega los datos extraidos con stratos a lo largo de la cuarta dimension (dias)
  }
  print("Datos de estratosfera extraidos")
  
  print('Combinando datos de troposfera y estratosfera')
  DATA[,,,,v,"2002"] = abind(tropos, stratos, along = 3) # concatena los datos de troposfera y estratosfera a lo largo de la dimension 3 (niveles de presion)
  print('Datos del 2002 extraidos')


##### DATOS 2019
  data2019 = NULL

  print('Estrayendo datos de 2019')
  for (m in months) { # Para cada mes
    nc = nc_open( # Abre el archivo nc
      here(DIR$ERA5_II, paste(variables[[v]][2], m, "2019_day.nc", sep = "_"))
    )
    data_extracted = ncvar_get( # Extrae los dato
      nc, # Archivo nc desde donde se van a extraer los datos
      varid = names(nc$var)[2], # Indica la variable a extraer
      start = c(
        1,
        which(ncvar_get(nc, "latitude") == lats[1]),
        1,1
        ), # Posicion de la primera longitud, latitud y nivel de presion
      count = c(-1,ny,8,-1) # Indica la cantidad de valores a leer a lo largo de cada dimension. Extrae de 1000 a 100 hPa.
    )
    nc_close(nc) # Cierra el archivo nc
    rm(nc) # Remueve la variable nc del espacio de trabajo
    
    data2019 = abind(data2019, data_extracted, along = 4)
   
  }
  
  DATA[,,,,v,"2019"] = data2019# Pega los datos extraidos en DATA a lo largo de la 4ta dimension (dias)
  print(paste('Datos de ', vars[v], ' extraidos.'))

}

rm(stratos,tropos,data2019,data_extracted)

########## Analasis Espectral #################################################
print("Analisis espectral...")
waves = c("p", "s"); nw = length(waves)

FOURIER = array(
  double(nx*ny*nz*nt*nv*na*nw), c(nx,ny,nz,nt,nv,na,nw),
  dimnames = list(lons = lons, lats = lats, plvl = plvl, dates = dates, vars = vars, years = years, waves = waves)
)

FOURIER[,,,,,,"p"] = apply(DATA, 2:6, FilterWave, k = 0:1)
FOURIER[,,,,,,"s"] = apply(DATA, 2:6, FilterWave, k = c(0,2))

########## Componente Vertical del Flujo de Eliassen Palm #####################
### Variables #################################
phi = lats*pi/180 # Coordenadas latitudinales en radianes
psr = plvl*100 # Presion en pascales

### Constantes ################################
r = 6.37122e06 # [m] Radio terrestre medio
omega = 7.2921e-5 # Rotacion Terrestre
ps = 100000 # [Pa] Presion de referencia estandarizada
R = 287 # [J.K-1.kg-1] Constante de los gases ideales
Cp = 1003.5 # [J.K-1.kg-1] Capacidad Calorifica del aire seco, al nivel del mar y a 273.15K
kpp = R/Cp # [adimensional] 

### PARAMETROS RELEVANTES #####################
f = 2*omega*sin(phi) # Parametro de Coriolis
acphi = r*cos(phi) # r por el coseno de la latitud
factor = array(f*acphi, c(ny,nz,nt,na,nw))

### TEMPERATURA POTENCIAL #####################
print("Calculando FEP...")
# Calculo Temperatura potencial
TP = FOURIER[,,,,"t",,]*array(rep((psr/ps)^-kpp, each = nx*ny), c(nx,ny,nz,nt,na,nw))
# Aplico media zonal
TP = apply(TP , 2:6, mean) 

### DERIVADA DE TP ############################
# Derivada de la media zonal de TP respecto a p
TP_p = array( # Array donde voy a guardar la derivada
  double(ny*nz*nt*na*nw), c(ny,nz,nt,na,nw), dimnames = list(lats = lats, plvl = plvl, dates = dates,years = years, waves = waves)
) 
for (z in 2:(nz-1)) { # Diferencias finitas centradas en dominio interior:
  TP_p[,z,,,] = (TP[,z+1,,,]-TP[,z-1,,,])/(psr[z+1]-psr[z-1])
}
TP_p[,1,,,] = (TP[,2,,,]-TP[,1,,,])/(psr[2]-psr[1])
TP_p[,nz,,,] = (TP[,nz,,,]-TP[,nz-1,,,])/(psr[nz]-psr[nz-1])

### Flujos de Calor ###########################
HF = array(
  double(ny*nz*nt*na*nw), c(ny,nz,nt,na,nw),
  dimnames = list(lats = lats, plvl = plvl, dates = dates,years = years, waves = waves)
)

for (w in 1:nw) {
  for (a in 1:na) {
    for (t in 1:nt) {
      for (z in 1:nz) {
        for (y in 1:ny)
          HF[y,z,t,a,w] = cov(FOURIER[,y,z,t,"t",a,w],FOURIER[,y,z,t,"v",a,w])
      }
    }
  }
}

### Componente Vertical del FEP ################
# Calculo del flujo:
Fp = factor*HF/TP_p

# # Escalado
# scale = 1/sqrt(1000/plvl)
# Fp = Fp*array(rep(scale, each = ny), c(ny,nz,nt,na,nw))

# Promedio latitud
Fp = apply(Fp, 2:5, mean)*1e-6



########## Elementos para graficar ###########################################3
print("Graficando...")
# mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
# mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
mini = c(213, 222, 232, 244, 253, 263, 273)
mname = c("1-Ago", "10", "20", "1-Sep", "10", "20", "30")

load("inputs/heatflux_mindays.RData")

neg = seq(-12,-1.5,1.5); nn = length(neg)
pos = seq(0,3,1.5); np = length(pos)
brks = c(neg,pos)
labs = as.character(brks)
labs[which(even(1:(nn+np)))] = " "

colors = c(
  sequential_hcl(nn, palette = "Blues", rev = FALSE), # YlGnBu / Blues /
  sequential_hcl(np-1, palette = "Reds", rev = TRUE) # YlOrRd / Reds / OrRd
)

#plot_title = list(p = "Ondas planetarias [k = 1 a 4]", s = "Ondas sinópticas [k = 5 a 12]")
plot_title = array(
  c("Onda-1", "Onda-2", "Onda-1", "Onda-2"),
  c(2,2), dimnames = list(c("p","s"), c("2002","2019"))
)
hfmin = list("2002" = mins2002, "2019" = mins2019)

tm.plot_margin = list(
  "p" = list(
    "2002" = margin(t = 0, r = 0.1, b = 0.01, l = 0, unit = "cm"), 
    "2019" = margin(t = 0.01, r = 0.1, b = 0, l = 0, unit = "cm")
  ),
  "s" = list(
    "2002" = margin(t = 0, r = 0, b = 0.01, l = 0.1, unit = "cm"), 
    "2019" = margin(t = 0.01, r = 0, b = 0, l = 0.1, unit = "cm")
  )
)
tm.plot_title = list("2002" = element_text(size = 7, vjust = -0.5, margin = margin(c(0,0,0,0), unit = "cm")), "2019" = element_blank())
tm.axis_text_x_bottom = list("2002" = element_text(), "2019" = element_blank())
tm.axis_ticks_x_bottom = list("2002" = element_line(), "2019" = element_blank())
tm.axis_ticks_x_top = list("2002" = element_blank(), "2019" = element_line())
tm.axis_text_y_left = list("p" = element_text(), "s" = element_blank())
tm.axis_text_y_right = list("p" = element_blank(), "s" = element_text())
tm.axis_title_y_left = list("p" = element_text(margin = margin(r = -0.01, unit = "cm")), "s" = element_blank())
tm.axis_title_y_right = list("p" = element_blank(), "s" = element_text(margin = margin(l = -0.01, unit = "cm")))
tm.axis_ticks_y_left = list("p" = element_line(), "s" = element_blank())
tm.axis_ticks_y_right = list("p" = element_blank(), "s" = element_line())

######### Grafico alternativa #################################################
GG = vector("list", nw*na)

MC = matrix(
  c("2002", "p",
    "2002", "s",
    "2019", "p",
    "2019", "s"),
  byrow = TRUE, ncol = 2, nrow = 4
)

for (i in 1:4) {
  a = MC[i,1]
  w = MC[i,2]
  
  df = data.frame(
    days = rep(days, nz),
    plvl = rep(plvl, each = nt),
    var = array(t(Fp[,,a,w]))
  )
  
  GG[[i]] = ggplot(df) +
    aes(x = days, y = plvl, z = var) +
    geom_contour_fill(breaks = brks) +
    scale_fill_gradientn(
      colors = colors,
      breaks = brks,
      labels = labs,
      limits = range(brks)
    ) +
    guides(
      fill = guide_colorstrip(
        title = expression("F"[p]~ "[10"^{6}~"Pa.m"^{2}~".s"^{-2}~"]    "), ticks = TRUE, ticks.color = "gray30", barwidth = 10, barheight = 0.3,
        label.vjust = 2
      )
    ) +
    geom_contour2(breaks = 0, size = 0.1, color = "gray40", alpha = 0.6) +
    scale_x_continuous(
      name = NULL, limits = range(213:273), breaks = mini, labels = mname, expand = c(0,0), sec.axis = dup_axis(name = NULL),
    ) +
    scale_y_level(
      name = "hPa", breaks = c(1000,850,700,500,400,300,200,100), sec.axis = dup_axis(name = "hPa")
    ) +
    geom_vline(xintercept = hfmin[[a]], color = "blue", linetype = "dashed", size = 0.1) +
    geom_hline(yintercept = c(500,200), color = "gray50", linetype = "dashed", size = 0.05) +
    annotate("text", x = 216, y = 900, label = a, color = "black", size = 1.75) +
    annotate("text", x = 217, y = 110, label = plot_title[w,a], color = "black", size = 2) +
    #labs(title = plot_title[[w]]) +
    theme_light() +
    theme(
      #plot.title = tm.plot_title[[a]],
      plot.title = element_blank(),
      plot.margin = tm.plot_margin[[w]][[a]],
      legend.title = element_text(size = 5, face = "bold"),
      legend.text = element_text(size = 4),
      legend.position = "bottom",
      legend.margin = margin(c(0,0,0,0), unit = "cm"),
      legend.box.spacing = unit(0, units = "cm"),
      axis.title = element_text(size = 5, face = "bold", color = "gray50"),
      axis.text.x.bottom = tm.axis_text_x_bottom[[a]],
      axis.text.x.top = element_blank(),
      axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[a]],
      axis.ticks.x.top = tm.axis_ticks_x_top[[a]],
      axis.title.y.left = tm.axis_title_y_left[[w]],
      axis.title.y.right = tm.axis_title_y_right[[w]],
      axis.text.y.left = tm.axis_text_y_left[[w]],
      axis.text.y.right = tm.axis_text_y_right[[w]],
      axis.ticks.y.left = tm.axis_ticks_y_left[[w]],
      axis.ticks.y.right = tm.axis_ticks_y_right[[w]],
      axis.text = element_text(size = 5)
    )
}

figure = ggarrange(plotlist = GG, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

print("Exportando...")
ggsave(
  here("outputs/ftropos/FEPp_vertical_tropos_v3.png"),
  figure, device = "png", width = 16, height = 9, units = "cm", dpi = 200
)
print("Hecho.")