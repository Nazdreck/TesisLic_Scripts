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
lats = seq(-45,-75,-2.5); ny = length(lats)
plvlt = c(500,400,300); nzt = length(plvlt)
plvls = c(200,100,70); nzs = length(plvls)
plvl = c(400,100); nz = length(plvl)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("01","02","03","04","05","06","07","08","09","10","11","12") 
vars = c("t","v"); nv = length(vars)
years = 1979:2019; na = length(years)

DATA = array(
  double(nx*ny*(nzs+nzt)*nt*na*nv), c(nx,ny,(nzs+nzt),nt,na,nv),
  dimnames = list(lons = lons, lats = lats, plvl = c(plvlt,plvls), days = days, years = years, vars = vars)
)

########## Archivos de Datos ##################################################
##### Lista de nombres de variables 
# Primera columna es para ERA5_I. Segunda columna es para era "ERA5_II"
variables = list( 
  Temp = c("temperature","t"), 
  Vwind = c("v_component_of_wind", "v")
)
########## Extraccion de datos ################################################
##### DATOS 1979-2018
### Troposfera
partial = NULL
print('Extrayendo datos 1979-2018')
for (v in 1:nv) {
  print(paste0('Extrayendo datos de ', vars[v]))
  nc_tropos = nc_open( # Abre el archivo nc que contiene los datos de la troposfera
    here(DIR$ERA5_I, paste0("era5.", variables[[v]][1], ".day.mean.nc"))
  )
  for (a in as.character(years[1:(na-1)])) {
    print(paste0('Extrayendo año ', a))
    if (!is.leapyear(as.numeric(a))) {
      ext = paste0(a, "-01-01")
      
      tropos = ncvar_get( # Extrae los datos de la troposfera
        nc_tropos, # Archivo nc desde donde se van a extraer los datos
        varid = names(nc_tropos$var)[2], # Indica la variable a extraer
        start = c( # Indica donde comenzar a leer
          1,
          which(ncvar_get(nc_tropos, "lat") == lats[1]),
          which(ncvar_get(nc_tropos, "lev") == plvlt[1]), # Posicion de la primera longitud, latitud y nivel de presion
          which(as.character(as.Date(ncvar_get(nc_tropos, "time")/24, origin = "1900-01-01")) == ext) # Posicion de la primera fecha
        ),
        count = c( # Indica la cantidad de valores a leer a lo largo de cada dimension
          -1,ny,nzt, # Lee todos los datos en las dimensiones x,y,z
          nt # Lee un año entero de datos, cambiar por 366 si el año es biciesto
        )
      )
      
      DATA[,,as.character(plvlt),,a,v] = tropos
    } else {
      ext = paste0(a, "-01-01")
      
      tropos1 = ncvar_get( # Extrae los datos de la troposfera
        nc_tropos, # Archivo nc desde donde se van a extraer los datos
        varid = names(nc_tropos$var)[2], # Indica la variable a extraer
        start = c( # Indica donde comenzar a leer
          1,
          which(ncvar_get(nc_tropos, "lat") == lats[1]),
          which(ncvar_get(nc_tropos, "lev") == plvlt[1]), # Posicion de la primera longitud, latitud y nivel de presion
          which(as.character(as.Date(ncvar_get(nc_tropos, "time")/24, origin = "1900-01-01")) == ext) # Posicion de la primera fecha
        ),
        count = c( # Indica la cantidad de valores a leer a lo largo de cada dimension
          -1,ny,nzt, # Lee todos los datos en las dimensiones x,y,z
          59 # Lee un año entero de datos, cambiar por 366 si el año es biciesto
        )
      )
      
      ext = paste0(a, "-03-01") 
      tropos2 = ncvar_get( # Extrae los datos de la troposfera
        nc_tropos, # Archivo nc desde donde se van a extraer los datos
        varid = names(nc_tropos$var)[2], # Indica la variable a extraer
        start = c( # Indica donde comenzar a leer
          1,
          which(ncvar_get(nc_tropos, "lat") == lats[1]),
          which(ncvar_get(nc_tropos, "lev") == plvlt[1]), # Posicion de la primera longitud, latitud y nivel de presion
          which(as.character(as.Date(ncvar_get(nc_tropos, "time")/24, origin = "1900-01-01")) == ext) # Posicion de la primera fecha
        ),
        count = c( # Indica la cantidad de valores a leer a lo largo de cada dimension
          -1,ny,nzt, # Lee todos los datos en las dimensiones x,y,z
          nt-59 # Lee un año entero de datos, cambiar por 366 si el año es biciesto
        )
      )
      DATA[,,as.character(plvlt),,a,v] = abind(tropos1,tropos2, along = 4)
    }
    
### Estratosfera

    for (m in months) {
      if (m == "02") s = 28 else s = -1
      nc = nc_open(
        here(DIR$ERA5_II, paste(vars[v],m,a, "day.nc", sep = "_"))
      )
      data_extracted = ncvar_get(
        nc,
        varid = names(nc$var)[2],
        start = c(
          1,
          which(ncvar_get(nc, "latitude") == lats[1]),
          which(ncvar_get(nc, "level") == plvls[1]),
          1
        ),
        count = c(
          -1,ny,nzs,s
        )
      )
      nc_close(nc); rm(nc)
      
      partial = abind(partial, data_extracted, along = 4)
    }
    DATA[,,as.character(plvls),,a,v] = partial
    partial = NULL
  }
}

nc_close(nc_tropos); rm(nc_tropos)

##### DATOS 2019
partial = NULL

print('Extrayendo datos 2019')
for (v in vars){
  for (p in c(850,200)) {
    for (m in months) {
      if (m == "02") s = 28 else s = -1
      nc = nc_open(
        here(DIR$ERA5_II, paste(v,m,"2019", "day.nc", sep = "_"))
      )
      data_extracted = ncvar_get(
        nc,
        varid = names(nc$var)[2],
        start = c(
          1,
          which(ncvar_get(nc, "latitude") == lats[1]),
          which(ncvar_get(nc, "level") == p),
          1
        ),
        count = c(
          -1,ny,3,s
        )
      )
      nc_close(nc); rm(nc)
      
      partial = abind(partial, data_extracted, along = 4)
    }
    
    if (p == 850) pp = plvlt else pp = plvls
    DATA[,,as.character(pp),,"2019",v] = partial
    partial = NULL
  }
}

rm(data_extracted, tropos, tropos1, tropos2, partial)
# ########## Flujos de Calor ####################################################
# HF = array(
#   double(ny*nz*nt*na), c(ny,nz,nt,na),
#   dimnames = list(lats = lats, plvl = plvl, days = days, years = years)
# )
# 
# ZZ = matrix(
#   c(2,400,
#     5,100),
#   byrow = TRUE, ncol = 2, nrow = 2
# )
# 
# for (a in 1:na) {
#   for (t in 1:nt) {
#     for (z in 1:nz) {
#       z1 = ZZ[z,1]
#       for (y in 1:ny) {
#         HF[y,z,t,a] = cov(DATA[,y,z1,t,a,"t"],DATA[,y,z1,t,a,"v"])
#       }
#     }
#   }
# }
# 
# ########## Perfiles de Theta ##################################################
# # Calculo Temperatura potencial
# theta = DATA[,,,,,"t"]*array(rep((c(plvlt,plvls)/1000)^-0.286, each = nx*ny), c(nx,ny,(nzs+nzt),nt,na))
# # Aplico media zonal
# zm_theta = apply(theta, 2:5, mean); rm(theta)
# 
# # Derivada de la media zonal de TP respecto a p
# logp = log(c(plvlt,plvls)) # Coordenadas log-P
# THETA_p = array( # Array donde voy a guardar la derivada
#   double(ny*nz*nt*na), c(ny,nz,nt,na), dimnames = list(lats = lats, plvl = plvl, dates = dates,years = years)
# ) 
# 
# 
# 
# for (i in 1:2) {
#   z1 = ZZ[i,1]
#   z2 = as.character(ZZ[i,2])
#   THETA_p[,z2,,] = ((zm_theta[,z1+1,,]-zm_theta[,z1-1,,])/(logp[z1+1]-logp[z1-1]))/(100*c(plvlt,plvls)[z1]) 
# }
# 
# ########## Componente Vertical del Flujo de Eliassen Palm #####################
# # Constantes
# r = 6.37122e06 # Radio terrestre
# omega = 7.2921e-5 # Rotacion Terrestre
# phi = lats*pi/180 # Coordenadas latitudinales en radianes
# f = 2*omega*sin(phi) # Parametro de Coriolis
# acphi = r*cos(phi) # r por el coseno de la latitud
# factor = array(f*acphi, c(ny,nz,nt,na))
# 
# # Calculo del flujo:
# Fp = factor*HF/THETA_p
# 
# # Escalado
# scale = 1/sqrt(plvl)
# Fp = Fp*array(rep(scale, each = ny), c(ny,nz,nt,na))
# 
# # Promedio latitud
# Fp = apply(Fp, 2:4, mean)/1000

########## Componente Vertical del Flujo de Eliassen Palm #####################
### Iterador ##################################
ZZ = matrix(
  c(2,400,
    5,100),
  byrow = TRUE, ncol = 2, nrow = 2
)

### Variables #################################
phi = lats*pi/180 # Coordenadas latitudinales en radianes
psr = c(plvlt,plvls)*100 # Presion en pascales

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
factor = array(f*acphi, c(ny,nz,nt,na))

### TEMPERATURA POTENCIAL #####################
print("Calculando FEP...")
# Calculo Temperatura potencial
TP = DATA[,,,,,"t"]*array(rep((psr/ps)^-kpp, each = nx*ny), c(nx,ny,nzt+nzs,nt,na))
# Aplico media zonal
TP = apply(TP , 2:5, mean) 

### DERIVADA DE TP ############################
# Derivada de la media zonal de TP respecto a p
TP_p = array( # Array donde voy a guardar la derivada
  double(ny*nz*nt*na), c(ny,nz,nt,na), dimnames = list(lats = lats, plvl = plvl, dates = dates,years = years)
) 
for (i in 1:2) { # Diferencias finitas centradas en dominio interior:
  z1 = ZZ[i,1]
  z2 = as.character(ZZ[i,2])
  TP_p[,z2,,] = (TP[,z1+1,,]-TP[,z1-1,,])/(psr[z1+1]-psr[z1-1])
}

### Flujos de Calor ###########################
HF = array(
  double(ny*nz*nt*na), c(ny,nz,nt,na),
  dimnames = list(lats = lats, plvl = plvl, days = days, years = years)
)
for (a in 1:na) {
  for (t in 1:nt) {
    for (z in 1:nz) {
      z1 = ZZ[z,1]
      for (y in 1:ny) {
        HF[y,z,t,a] = cov(DATA[,y,z1,t,a,"t"],DATA[,y,z1,t,a,"v"])
      }
    }
  }
}

### Componente Vertical del FEP ################
# Calculo del flujo:
Fp = factor*HF/TP_p

# Promedio latitud
Fp = apply(Fp, 2:4, mean)*1e-6

########## Percentiles ########################################################
pers = c(10,30,70,90); np = length(pers)

PER = array(double(np*nt*nz), c(nt,np,nz), dimnames = list(days = days,percentile = pers, plvl = plvl))

for (z in 1:2) {
  for (d in days) {
    for (p in pers) {
      x1 = which(Percentile(Fp[z,d,]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
      x2 = min(Fp[z,d,][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
      PER[d,as.character(p),z] = x2 # Guarda el valor minimo encontrado
    }
  }
}

########## Elementos para graficar ############################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")

brks = list(seq(-24,8,4), seq(-4.5,1.5,0.75))
leg = list("none",c(0.015, 0.045))
tit = list(expression("[400 hPa]"), expression("[100 hPa]"))
#axsx = list(element_blank(), element_text())
tm.axis_text_x_bottom = list("top" = element_text(), "bot" = element_blank())
tm.axis_ticks_x_bottom = list("top" = element_line(), "bot" = element_blank())
tm.axis_ticks_x_top = list("top" = element_blank(), "bot" = element_line())
tm.plot_margin = list(
  "top" = margin(t = 0.12, r = 0.2, b = 0, l = 0.27, unit = "cm"),
  "bot" = margin(t = 0.12, r = 0.2, b = 0, l = 0.1, unit = "cm")
)
########## Graficos ###########################################################
GG_amp = vector("list", 2)

for (k in 1:2) {
  
  df_lines = data.frame(
    days = rep(days, 9),
    vars = c(
      p10 = PER[,"10",k],
      p90 = PER[,"90",k],
      p30 = PER[,"30",k],
      p70 = PER[,"70",k],
      max = apply(Fp[k,,-c(24,41)], 1, max), 
      min = apply(Fp[k,,-c(24,41)], 1, min),
      clim = apply(Fp[k,,-c(24,41)], 1, mean),
      d2002 = Fp[k,,"2002"],
      d2019 = Fp[k,,"2019"]
    ),
    group = rep(LETTERS[1:9], each = nt),
    cols = c(
      rep("10-90%", each = nt*2), rep("30-70%", each = nt*2), rep("Máx/mín.", each = nt*2),
      rep("Clim.", each = nt), rep("2002", each = nt), rep("2019", each = nt)
    )
  )
  
  df_lines$cols = factor(df_lines$cols, levels = c("10-90%", "30-70%", "Máx/mín.", "Clim.", "2002", "2019"))
  
  df_ribbons = data.frame(
    days = rep(days, 2),
    pmin = c(p10 = PER[,"10",k], p30 = PER[,"30",k]),
    pmax = c(p90 = PER[,"90",k], p70 = PER[,"70",k]),
    group = rep(LETTERS[1:2], each = nt)
  )
  
  GG_amp[[k]] = ggplot() +
    geom_ribbon(
      data = df_ribbons, aes(x = days, ymin = pmin, ymax = pmax, fill = group), 
      alpha = 0.5, show.legend = FALSE
    ) +
    scale_fill_manual(values = c("A"= "gray80", "B" = "gray40")) +
    geom_line(
      data = df_lines, aes(x = days, y = vars, group = group, color = cols, size = cols)
    ) +
    scale_color_manual(
      name = expression("F"[p]~ "[10"^{6}~"Pa.m"^{2}~".s"^{-2}~"]"),
      values = c(
        "Máx/mín." = "darkgreen", "Clim." = "black", 
        "2019" = "blue", "2002"= "red", 
        "10-90%"= "gray70", "30-70%" = "gray30"
      )
    ) +
    scale_size_manual(
      values = c(
        "Máx/mín." = 0.15, "Clim." = 0.3, 
        "2019" = 0.3, "2002"= 0.3, 
        "10-90%"= 0.1, "30-70%" = 0.1
      )
    ) +
    guides(size = "none") +
    scale_x_continuous(name = NULL, limits = range(1:365), breaks = mini, labels = mname, expand = c(0,0), sec.axis = dup_axis(name = NULL))+
    scale_y_continuous(
      name = NULL, breaks = brks[[k]], limits = range(brks[[k]]),
      expand = c(0,0)
    ) +
    labs(
      title = tit[[k]]
    ) +
    theme_light() +
    theme(
      plot.margin = tm.plot_margin[[k]],
      plot.title = element_text(face = "bold", vjust = -9, hjust = 0.02, size = 7, margin = margin(t = -10)),
      legend.key.size = unit(0.275, "cm"),
      legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
      legend.margin = margin(-2, 3, 3, 3),
      legend.key = element_rect(fill = alpha("orange", 0.05)),
      legend.title = element_text(size = 5, face = "bold", vjust = -1),
      legend.position = leg[[k]],
      legend.justification = c("left", "bottom"),
      legend.text = element_text(size = 6),
      #axis.text.x = axsx[[k]],
      axis.text.x.bottom = tm.axis_text_x_bottom[[k]],
      axis.text.x.top = element_blank(),
      axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[k]],
      axis.ticks.x.top = tm.axis_ticks_x_top[[k]],
      axis.text = element_text(size = 6)
    )
}

fig_amp = ggarrange(
  plotlist = list(GG_amp[[1]],GG_amp[[2]]), ncol = 1, nrow = 2#, align = "hv"
)

ggsave(
  here("outputs/ftropos", "FEPp_series5.png"),
  fig_amp, device = "png", width = 12, height = 9, units = "cm", dpi = 200
)


###############################################################################
########## Promedio 45 dias ###################################################
Fp_45dm = array(
  double(nt*na), c(nt,na), dimnames = list(days = days, years = years)
)

for (a in 1:na) {
  for (t in 46:nt) {
    Fp_45dm[t,a] = mean(Fp["100",(t-45):(t-1),a])
  }
}
Fp_45dm[1:45,] = NA

########## Percentiles ########################################################
pers = c(10,30,70,90); np = length(pers)

PER = array(double(np*nt), c(nt,np), dimnames = list(days = days,percentile = pers))


for (d in days[46:nt]) {
  for (p in pers) {
    x1 = which(Percentile(Fp_45dm[d,]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
    x2 = min(Fp_45dm[d,][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
    PER[d,as.character(p)] = x2 # Guarda el valor minimo encontrado
  }
}

PER[1:45,] = NA
######### Grafico #############################################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
brks = seq(-1.2,0,0.2)


df_lines = data.frame(
  days = rep(days, 9),
  vars = c(
    p10 = PER[,"10"],
    p90 = PER[,"90"],
    p30 = PER[,"30"],
    p70 = PER[,"70"],
    max = apply(Fp_45dm[,-c(24,41)], 1, max),
    min = apply(Fp_45dm[,-c(24,41)], 1, min),
    clim = apply(Fp_45dm[,-c(24,41)], 1, mean),
    d2002 = Fp_45dm[,"2002"],
    d2019 = Fp_45dm[,"2019"]
  ),
  group = rep(LETTERS[1:9], each = nt),
  cols = c(
    rep("10-90%", each = nt*2), rep("30-70%", each = nt*2), rep("Máx/mín.", each = nt*2),
    rep("Clim.", each = nt), rep("2002", each = nt), rep("2019", each = nt)
  )
)

df_lines$cols = factor(df_lines$cols, levels = c("10-90%", "30-70%", "Máx/mín.", "Clim.", "2002", "2019"))

df_ribbons = data.frame(
  days = rep(days, 2),
  pmin = c(p10 = PER[,"10"], p30 = PER[,"30"]),
  pmax = c(p90 = PER[,"90"], p70 = PER[,"70"]),
  group = rep(LETTERS[1:2], each = nt)
)

GG = ggplot() +
  geom_ribbon(
    data = df_ribbons, aes(x = days, ymin = pmin, ymax = pmax, fill = group), 
    alpha = 0.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = c("A"= "gray80", "B" = "gray40")) +
  geom_line(
    data = df_lines, aes(x = days, y = vars, group = group, color = cols, size = cols)
  ) +
  scale_color_manual(
    name = expression("F"[p]~ "[10"^{6}~"Pa.m"^{2}~".s"^{-2}~"]"),
    values = c(
      "Máx/mín." = "darkgreen", "Clim." = "black", 
      "2019" = "blue", "2002"= "red", 
      "10-90%"= "gray70", "30-70%" = "gray30"
    )
  ) +
  scale_size_manual(
    values = c(
      "Máx/mín." = 0.15, "Clim." = 0.3, 
      "2019" = 0.3, "2002"= 0.3, 
      "10-90%"= 0.1, "30-70%" = 0.1
    )
  ) +
  guides(size = "none") +
  scale_x_continuous(name = NULL, limits = range(46:365), breaks = mini, labels = mname, expand = c(0,0))+
  scale_y_continuous(
    name = NULL, breaks = brks, limits = range(brks),
    expand = c(0,0)
  ) +
  labs(
    title = expression("F"[p]~ "(45-75°S, 100 hPa, 45 dias previos)")
  ) +
  theme_light() +
  theme(
    legend.key.size = unit(0.275, "cm"),
    legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
    legend.margin = margin(-2, 3, 3, 3),
    legend.key = element_rect(fill = alpha("orange", 0.05)),
    legend.title = element_text(size = 5, face = "bold", vjust = -1.6),
    legend.position = c(0.01,0.025),
    legend.justification = c("left", "bottom"),
    plot.title = element_text(face = "bold", vjust = -7.5, hjust = 0.02, size = 7, margin = margin(t = -10)),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
  )


ggsave(
  here("outputs/ftropos", "Fp_45dm_series4.png"),
  GG, device = "png", width = 12, height = 5.0625, units = "cm", dpi = 200
)



























