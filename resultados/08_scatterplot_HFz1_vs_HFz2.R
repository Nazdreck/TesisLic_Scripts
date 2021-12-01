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

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(-45,-75,-2.5); ny = length(lats)
#plvl = c(100,70,30,10); nz = length(plvl)
plvl = c(100,10); nz = length(plvl)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("01","02","03","04","05","06","07","08","09","10","11","12") 
vars = c("t","u","v"); nv = length(vars)
years = 1979:2019; na = length(years)

DATA = array(
  double(nx*ny*nz*nt*na*nv), c(nx,ny,nz,nt,na,nv),
  dimnames = list(lons = lons, lats = lats, plvl = plvl, days = days, years = years, vars = vars)
)

########## Directories ########################################################
##### Guarda en una lista las rutas de los datos
DIR = list(
  ALLDATA = "../../../datos", # Todos los datos de ERA5 en Server Vegeta
  ERA5_I = "../../../datos/ERA5/day", # Datos diarios para los 12 meses del año de T,U,V,Z,PV desde 1000 a 300 hpa entre 1979-2018.
  ERA5_II = "../../../datos/ERA5_aux/nuevos/daily" # Datos diarios para los 12 meses del año de T,U,V,Z,PV desde 200 a 1 hpa entre 1979-2018 y desde 1000 a 1 hpa para 2019
)

########## Extrayendo datos ###################################################
partial = NULL

for (v in vars){
  print(paste0("Extrayendo datos de ", v))
  for (a in as.character(years)) {
    print(paste0("Extrayendo año ", a))
    for (p in plvl) {
      for (m in months) {
        if (m == "02") s = 28 else s = -1
        nc = nc_open(
          here(DIR$ERA5_II, paste(v,m,a, "day.nc", sep = "_"))
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
            -1,ny,1,s
          )
        )
        nc_close(nc); rm(nc)
        
        partial = abind(partial, data_extracted, along = 3)
      }
      
      DATA[,,as.character(p),,a,v] = partial
      partial = NULL
    }
  }
}

########## ONDA 1 VS ONDA 2 ###################################################
########## Analisis Espectral #################################################
waves = c("k1", "k2"); nw = length(waves)
dd = 121:273; ndd = length(dd)

FOURIER = array(
  double(nx*ny*nz*ndd*(nv-1)*na*nw), c(nx,ny,nz,ndd,na,nv-1,nw),
  dimnames = list(lons = lons, lats = lats, plvl = plvl, days = dd, years = years, vars = vars[c(1,3)], waves = waves)
)

FOURIER[,,,,,,"k1"] = apply(DATA[,,,dd,,vars[c(1,3)]], 2:6, FilterWave, k = 1)
FOURIER[,,,,,,"k2"] = apply(DATA[,,,dd,,vars[c(1,3)]], 2:6, FilterWave, k = 2)

########## Flujos de Calor ####################################################
HF = array(
  double(ny*nz*ndd*na*nw), c(ny,nz,ndd,na,nw),
  dimnames = list(lats = lats, plvl = plvl, days = dd,years = years, waves = waves)
)

for (w in 1:nw) {
  for (a in 1:na) {
    for (t in 1:ndd) {
      for (z in 1:nz) {
        for (y in 1:ny) {
          HF[y,z,t,a,w] = cov(FOURIER[,y,z,t,a,"t",w],FOURIER[,y,z,t,a,"v",w])
        }
      }
    }
  }
}

HF_lm = apply(HF, 2:5, mean)


########## Elementos para graficar ############################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
brksx = list(
  "100" = seq(-60,30,30),
  # "70" = seq(-90,30,30),
  # "30" = seq(-200,100,100),
  "10" = seq(-400,100,100)
)
brksy = list(
  "100" = seq(-60,30,30),
  # "70" = seq(-60,40,20),
  # "30" = seq(-120,60,30),
  "10" = seq(-180,90,90)
)

leg = list(
  "100" = c(0.020, 0.015), 
  # "70" = "none","30" = "none",
  "10" = "none")

aa = as.character(c(1979:2001,2003:2018))

########## Grafico ############################################################
GG = vector("list", nz); names(GG) = plvl

for (p in as.character(plvl)) {
  df_gen = data.frame(
    k1 = array(HF_lm[p,,aa,"k1"]),
    k2 = array(HF_lm[p,,aa,"k2"])
  )
  
  df_ssw = data.frame(
    k1 = array(HF_lm[p,,c("2002","2019"),"k1"]),
    k2 = array(HF_lm[p,,c("2002","2019"),"k2"]),
    group = c(rep("2002",ndd), rep("2019",ndd))
  )
  
  GG[[p]] = ggplot() +
    geom_point(
      data = df_gen, aes(x = k1, y = k2),
      color = "gray40", alpha = 0.5, size = 0.4
    ) +
    geom_point(
      data = df_ssw, aes(x = k1, y = k2, group = group, color = group),
      size = 0.6, alpha = 0.8
    ) +
    scale_color_manual(
      name = "Leyenda",
      values = c(
        "2019" = "blue", "2002"= "red"
      )
    ) +
    geom_hline(yintercept = 0, size = 0.2, linetype = "dashed", color = "gray15") +
    geom_vline(xintercept = 0, size = 0.2, linetype = "dashed", color = "gray15") +
    scale_x_continuous(
      name = NULL, expand = c(0,0), breaks = brksx[[p]], limits = range(brksx[[p]])
    ) +
    scale_y_continuous(
      name = NULL, expand = c(0,0), breaks = brksy[[p]], limits = range(brksy[[p]])
    ) +
     labs(
      title = paste0(p, " hPa")
    ) +
    theme_light() +
    theme(
      plot.title = element_text(face = "bold", vjust = -7, hjust = 0.08, size = 6, margin = margin(t = -8)),
      axis.text = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.justification = c("left", "bottom"),
      legend.position = leg[[p]],
      legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
      legend.margin = margin(-4, 1, 1, 1),
      legend.key = element_rect(fill = alpha("orange", 0.05)),
      legend.key.size = unit(0.275, "cm")
    )
}

figure = ggarrange(
  plotlist = GG, nrow = 1, ncol = 2, common.legend = FALSE
)

figure = annotate_figure(
  figure,
  left = text_grob( expression("Transporte "~ bar("v'T'")~ " por onda-2 [K ms"^{-1}*"]"), size = 6, rot = 90),
  bottom = text_grob( expression("Transporte "~ bar("v'T'")~ " por onda-1 [K ms"^{-1}*"]"), size = 6)
)

ggsave(
  here("outputs/ftropos","HF_scatter_k1_k2_p_v2.png"),
  figure, device = "png", width = 9, height = 6.75, units = "cm", dpi = 200
)

######### Flujos de Calor vs viento zonal #####################################
########## Viento Zonal #######################################################
dd = 121:273; ndd = length(dd)
U_zm = apply(DATA[,,,dd,,"u"], 2:5, mean)

########## Transporte de Calor ################################################
HFt = array(
  double(ny*nz*ndd*na), c(ny,nz,ndd,na),
  dimnames = list(lats = lats, plvl = plvl, days = dd,years = years)
)


for (a in 1:na) {
  for (t in 1:ndd) {
    for (z in 1:nz) {
      for (y in 1:ny) {
        HFt[y,z,t,a] = cov(DATA[,y,z,t,a,"t"],DATA[,y,z,t,a,"v"])
      }
    }
  }
}


HFt_lm = apply(HFt, 2:4, mean)

########## Elementos para graficar ############################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
brksx = list(
  "100" = seq(10,40,10),
  "70" = seq(10,50,10),
  "30" = seq(-15,60,15),
  "10" = seq(-20,100,20)
)
brksy = list(
  "100" = seq(-45,15,15),
  "70" = seq(-45,15,15),
  "30" = seq(-75,15,15),
  "10" = seq(-150,25,25)
)

leg = list("100" = c(0.020, 0.015), "70" = "none","30" = "none","10" = "none")

aa = as.character(c(1979:2001,2003:2018))

######### Grafico #############################################################
GG = vector("list", nz); names(GG) = plvl
ll = "-45"
p = "100"

for (p in as.character(plvl)) {
  df_gen = data.frame(
    u = array(U_zm[ll,p,,aa]),
    hf = array(HFt_lm[p,,aa])
  )
  
  df_ssw = data.frame(
    u = array(U_zm[ll,p,,c("2002","2019")]),
    hf = array(HFt_lm[p,,c("2002","2019")]),
    group = c(rep("2002",ndd), rep("2019",ndd))
  )
  
  GG[[p]] = ggplot() +
    geom_point(
      data = df_gen, aes(x = u, y = hf),
      color = "gray40", alpha = 0.5, size = 0.4
    ) +
    geom_point(
      data = df_ssw, aes(x = u, y = hf, group = group, color = group),
      size = 0.6, alpha = 0.8
    ) +
    scale_color_manual(
      name = "Leyenda",
      values = c(
        "2019" = "blue", "2002"= "red"
      )
    ) +
    geom_hline(yintercept = 0, size = 0.2, linetype = "dashed", color = "gray15") +
    geom_vline(xintercept = 0, size = 0.2, linetype = "dashed", color = "gray15") +
    scale_x_continuous(
      name = NULL, expand = c(0,0), breaks = brksx[[p]], limits = range(brksx[[p]])
    ) +
    scale_y_continuous(
      name = NULL, expand = c(0,0), breaks = brksy[[p]], limits = range(brksy[[p]])
    ) +
    labs(
      title = paste0(p, " hPa")
    ) +
    theme_light() +
    theme(
      plot.title = element_text(face = "bold", vjust = -7, hjust = 0.08, size = 6, margin = margin(t = -8)),
      axis.text = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      legend.justification = c("left", "bottom"),
      legend.position = leg[[p]],
      legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
      legend.margin = margin(-4, 1, 1, 1),
      legend.key = element_rect(fill = alpha("orange", 0.05)),
      legend.key.size = unit(0.275, "cm")
    )
}

figure = ggarrange(
  plotlist = GG, nrow = 1, ncol = 4, common.legend = FALSE
)

figure = annotate_figure(
  figure,
  left = text_grob( expression("Transporte "~ bar("v'T'")~ " [K ms"^{-1}*"]"), size = 6, rot = 90),
  bottom = text_grob( bquote(bar("U")~" en "~.(-as.numeric(ll))~"°S [ms"^{-1}*"]"), size = 6)
)

ggsave(
  here("outputs/ftropos",paste0("scatter_HF_U_p_",-as.numeric(ll),"S.png")),
  figure, device = "png", width = 16, height = 6.75, units = "cm", dpi = 200
)













