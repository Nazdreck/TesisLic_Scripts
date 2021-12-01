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
lats = seq(-30,-90,-2.5); ny = length(lats)
plvl = c(100,70,30,10); nz = length(plvl)
days = 152:243; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("06","07","08") 
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

rm(data_extracted,partial)

########## Transporte de Calor ################################################
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


HFyd_mean = apply(HF, c(2,4), mean)
########## Viento Zonal #######################################################
U_zm = apply(DATA[,,,,,"u"], c(2,3,5), mean)

########## Elementos para graficar ############################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
aa = as.character(c(1979:2001,2003:2018))
p = "10"
lts = c("-45", "-70")


tm.plot_margin = list(
  "-45" = margin(t = 0.2, r = 0, b = 0, l = 0.2, unit = "cm"), 
  "-70" = margin(t = 0.2, r = 0.2, b = 0, l = 0, unit = "cm")
)
# margin = margin(r = 0.1, unit = "cm")
tm.axis_title_y_left = list("-45" = element_text(), "-70" = element_blank())
tm.axis_text_y_left = list("-45" = element_text(), "-70" = element_blank())
tm.axis_ticks_y_left = list("-45" = element_line(), "-70" = element_blank())

########## Graficos ###########################################################
GG = vector("list", 2); names(GG) = c("-45", "-70")

for (ll in lts) {
  df_gen = data.frame(
    u = array(U_zm[ll,p,aa]),
    hf = (10/100)*array(HFyd_mean["10",aa]/HFyd_mean["100",aa]),
    year = substr(aa,3,4)
  )
  
  df_ssw = data.frame(
    u = array(U_zm[ll,p,c("2002","2019")]),
    hf = (10/100)*array(HFyd_mean["10",c("2002", "2019")]/HFyd_mean["100",c("2002", "2019")]),
    year = substr(c("2002", "2019"),3,4),
    group = c("2002", "2019")
  )
  
  GG[[ll]] = ggplot() +
    geom_text(
      data = df_gen, aes(x = u, y = hf, label = year), size = 2
    ) +
    geom_text(
      data = df_ssw, aes(x = u, y = hf, group = group, color = group, label = year), size = 3
    ) +
    scale_color_manual(
      name = "Leyenda",
      values = c(
        "2019" = "blue", "2002"= "red"
      )
    ) +
    guides(
      color = "none"
    ) +
    scale_x_continuous(
      name = bquote(bar("u")~" en "~.(-as.numeric(ll))~"°S [ms"^{-1}*"]")#, expand = c(0,0), breaks = seq(0,35,5), limits = range(seq(0,25,5))
    ) +
    scale_y_continuous(
      expression("Cociente "~ bar("v'T'")), expand = c(0,0),breaks = seq(0.2,0.8,0.1), limits = range(seq(0.2,0.8,0.1))
    ) +
    theme_light() +
    theme(
      plot.margin = tm.plot_margin[[ll]],
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.title.y.left = tm.axis_title_y_left[[ll]],
      axis.text.y.left = tm.axis_text_y_left[[ll]],
      axis.ticks.y.left = tm.axis_ticks_y_left[[ll]]
    )
}

figure = ggarrange(
  plotlist = GG, ncol = 2, nrow = 1
)

ggsave(
  here("outputs/ftropos",paste0("scatter_year_HF_u3.png")),
  figure, device = "png", width = 16, height = 6.75, units = "cm", dpi = 200
)

########## Posicion del maximo del jet ######################################## 
JMP = array(double(nz*na), c(nz,na), dimnames = list(plvl = plvl, years = years))

for (a in 1:na) {
  JMP["100",a] = lats[which(U_zm[,"100",a] == max(U_zm[,"100",a]))]
  JMP["10",a] = lats[which(U_zm[,"10",a] == max(U_zm[,"10",a]))]
}

JMP_clim = apply(JMP, 1, mean)
clim
########## Elementos para graficar ############################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
aa = as.character(c(1979:2001,2003:2018))

pls = c("100","10")

tm.plot_margin = list(
  "100" = margin(t = 0.2, r = 0, b = 0, l = 0.2, unit = "cm"), 
  "10" = margin(t = 0.2, r = 0.2, b = 0, l = 0, unit = "cm")
)
# margin = margin(r = 0.1, unit = "cm")
tm.axis_title_y_left = list("100" = element_text(), "10" = element_blank())
tm.axis_text_y_left = list("100" = element_text(), "10" = element_blank())
tm.axis_ticks_y_left = list("100" = element_line(), "10" = element_blank())

########## Graficos ###########################################################
GG = vector("list", 2); names(GG) = pls

for (pp in pls) {
  df_gen = data.frame(
    mlats = array(JMP[pp,aa]),
    hf = (10/100)*array(HFyd_mean["10",aa]/HFyd_mean["100",aa]),
    year = substr(aa,3,4)
  )
  
  df_ssw = data.frame(
    mlats = array(JMP[pp,c("2002","2019")]),
    hf = (10/100)*array(HFyd_mean["10",c("2002", "2019")]/HFyd_mean["100",c("2002", "2019")]),
    year = substr(c("2002", "2019"),3,4),
    group = c("2002", "2019")
  )
  
  GG[[pp]] = ggplot() +
    geom_text(
      data = df_gen, aes(x = mlats, y = hf, label = year), size = 2
    ) +
    geom_text(
      data = df_ssw, aes(x = mlats, y = hf, group = group, color = group, label = year), size = 3
    ) +
    annotate("text", x= )
    scale_color_manual(
      name = "Leyenda",
      values = c(
        "2019" = "blue", "2002"= "red"
      )
    ) +
    guides(
      color = "none"
    ) +
    scale_x_continuous(
      name = bquote("Posicíon del máximo del jet en "~.(-as.numeric(pp))~"hPa")#, expand = c(0,0), breaks = seq(0,35,5), limits = range(seq(0,25,5))
    ) +
    scale_y_continuous(
      expression("Cociente "~ bar("v'T'")), expand = c(0,0),breaks = seq(0.2,0.8,0.1), limits = range(seq(0.2,0.8,0.1))
    ) +
    theme_light() +
    theme(
      plot.margin = tm.plot_margin[[pp]],
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.title.y.left = tm.axis_title_y_left[[pp]],
      axis.text.y.left = tm.axis_text_y_left[[pp]],
      axis.ticks.y.left = tm.axis_ticks_y_left[[pp]]
    )
}

figure = ggarrange(
  plotlist = GG, ncol = 2, nrow = 1
)

ggsave(
  here("outputs/ftropos",paste0("scatter_year_HF_jmp2.png")),
  figure, device = "png", width = 16, height = 6.75, units = "cm", dpi = 200
)













