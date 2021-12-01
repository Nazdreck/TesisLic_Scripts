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
require(maps)           # Draw geographical maps
require(mapproj)        # Map Projections
require(ggtext)
require(gtools)

########## Predefiniendo Variables ############################################
### Baja Resolucion
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(-12.5,-90,-2.5); ny = length(lats)
### Alta resolucion
lons_hr = seq(-180,179.75,0.25); nxhr = length(lons_hr)
lats_hr = seq(-12.5,-90,-0.25); nyhr = length(lats_hr)
###
stat = c("2002","2019","climatology"); ns = length(stat)
years = c("2002", "2019"); na = length(years)
days = 1:92; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "1900-08-01")), 6,10)
vars = c("t","gh"); nv = length(vars)
plvl = 10

DATA = array(
  double(nx*ny*nt*ns*nv), c(nx,ny,nt,ns,nv), 
  dimnames = list(lons = lons, lats = lats, dates = dates, stat = stat, vars = vars)
)

PV = array(
  double(nxhr*nyhr*nt*na), c(nxhr, nyhr, nt, na),
  dimnames = list(lons = lons_hr, lats = lats_hr, dates = dates, years = years)
)

########## Extrayendo datos ###################################################
print("Extrayendo datos...")
for (v in vars) {
  for (s in stat) {
    nc = nc_open(
      here("inputs/era5_plvls/", paste0(v, "_", s, ".nc"))
    )
    for (z in plvl){ 
      data_extracted = ncvar_get(
        nc,
        start = c(
          1,
          which(ncvar_get(nc, "lats") == lats[1]),
          which(ncvar_get(nc, "plvl") == plvl),
          which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "08-01")
        ),
        count = c(
          -1,ny,1,nt
        )
      )
      
      
      DATA[,,,s,v] = data_extracted
    }
    nc_close(nc); rm(nc)
  }
}
rm(data_extracted)

partial = NULL
for (a in as.character(years)) {
  for (m in c("08","09","10")) {
    nc = nc_open(
      here("inputs/era5_potvot/daily", paste0("pv_",m, "_", a,"_day.nc"))
    )
    
    partial = abind(partial, ncvar_get(
      nc, 
      varid = names(nc$var)[2],
      start = c(1,which(ncvar_get(nc, "latitude") == lats_hr[1]),1),
      count = c(-1,nyhr,-1)
    ), along = 3)
    nc_close(nc); rm(nc)
  }
  PV[,,,a] = partial
  partial = NULL
}
rm(partial)


print("Extraido.")
########## Acomodando los datos ###############################################
##### Extrayendo dias a graficar
tt = c(-31,-26,-21,-16,-11,-6,-4,-2,0,2,5,12); ntt = length(tt)
dd = list("2002" = dates[which(dates == "09-27") + tt], "2019" = dates[which(dates == "09-19") + tt])

extDAYs = array(
  double(nx*ny*ntt*na*nv), c(nx,ny,ntt,na,nv), 
  dimnames = list(lons = lons, lats = lats, dates = tt, years = years, vars = vars)
)
pv_extDAYs = array(
  double(nxhr*nyhr*ntt*na), c(nxhr,nyhr,ntt,na), 
  dimnames = list(lons = lons_hr, lats = lats_hr, dates = tt, years = years)
)
extDAYs[,,,"2002",] = DATA[,,dd[["2002"]],"2002",]
extDAYs[,,,"2019",] = DATA[,,dd[["2019"]],"2019",]

pv_extDAYs[,,,"2002"] = PV[,,dd[["2002"]], "2002"]
pv_extDAYs[,,,"2019"] = PV[,,dd[["2019"]], "2019"]

##### Acomodando los datos para proyeccion polar
### 2002 y 2019
d360 = extDAYs["180",,,,]
d360 = abind(d360, extDAYs[which(lons > 180) ,,,,], along = 1)
d360 = abind(d360, extDAYs[which(lons <= 180),,,,], along = 1)

pv_d360 = pv_extDAYs["-180",,,]
pv_d360 = abind(pv_extDAYs, pv_d360, along = 1)
# pv_d360 = abind(pv_d360, pv_extDAYs[which(lons_hr > 180) ,,,], along = 1)
# pv_d360 = abind(pv_d360, pv_extDAYs[which(lons_hr <= 180),,,], along = 1)

### climatologia
c360 = DATA["180",,,"climatology",]
c360 = abind(c360, DATA[which(lons > 180) ,,,"climatology",], along = 1)
c360 = abind(c360, DATA[which(lons <= 180),,,"climatology",], along = 1)

### nombrando arrays
lons = seq(-180,180,2.5); nx = length(lons)
dimnames(d360) = list(lons = lons, lats = lats, dates = tt, years = years, vars = vars)
dimnames(c360) = list(lons = lons, lats = lats, dates = dates, vars = vars)

lons_hr = seq(-180,180,0.25); nxhr = length(lons_hr)
dimnames(pv_d360) =  list(lons = lons_hr, lats = lats_hr, dates = tt, years = years)

##### cambio de unidades
d360[,,,,"gh"] = d360[,,,,"gh"]/9800
c360[,,,"gh"] = c360[,,,"gh"]/9800
pv_d360 = pv_d360/10^-6


########## Elementos para graficar ############################################
##### TITULOS
mon = list("2002" = substr(dd[["2002"]],1,2),"2019" = substr(dd[["2019"]],1,2))
for (a in years) {
  for (t in 1:ntt) {
    mon[[a]][t] = switch (as.numeric(mon[[a]][t])-7,"AGO", "SEP", "OCT")
  }
}

ttl = list(
  "2002" = paste(substr(dd[["2002"]],4,5), mon[["2002"]], "2002", sep = "-"),
  "2019" = paste(substr(dd[["2019"]],4,5), mon[["2019"]], "2019", sep = "-")
)

stl = c("-31d","-26d","-21d","-16d","-11d","-6d","-4d","-2d","0d","+2d","+5d","+12d")

##### QUIEBRES
brks = list(
  "t"  = seq(-20,60,10),
  "gh" = seq(26.5,32.5,0.5),
  "pv" = seq(-1380,320,20) # otra opcion: seq(-1300,200,100)
)
nb = list("t" =  length(brks$t),"gh" = length(brks$gh),"pv" = length(brks$pv))

##### LABELS
labs = vector("list", nv); names(labs) = vars
labs$t = as.character(brks$t); labs$t[which(odd(1:nb$t))] = " "
labs$gh = as.character(brks$gh); labs$gh[which(odd(1:nb$gh))] = " "
labs$pv = as.character(brks$pv)#; labs$pv[which(even(1:nb$pv))] = " "

for (b in 1:nb$pv) {
  if (!any(as.character(c(-1300,-1000,-700,-400,-100,300)) == labs$pv[b])) {
    labs$pv[b] = " "
  }
}

##### PV Colors
pv_cols = c(
  sequential_hcl(38, palette = "RdPu")[1:19],
  sequential_hcl(30, palette = "GnBu")[1:15],
  sequential_hcl(30, palette = "YlGn")[1:15],
  sequential_hcl(30, palette = "YlOrRd")[1:15],
  sequential_hcl(42, palette = "Turku", rev = TRUE)[22:42]#,
  #sequential_hcl(42, palette = "Lajolla", rev = TRUE)[7:12]
)

t_cols = c(
  sequential_hcl(10, palette = "Blues")[1:2],
  "#FFFFFF",
  sequential_hcl(6, palette = "PinkYl", rev = TRUE)
)

#####  MAPA
SH = map_data("world", ylim = c(-90,lats[1]))

##### MATRIZ iteradora

A = matrix(
  c(c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,7,8,9,10,9,10,11,12,11,12),
    rep(c(2002,2002,2019,2019), 6)),
  ncol = 2, nrow= 24, byrow = FALSE
)

#### anomalia de temperatura
anom_t = d360[,,,,"t"]
for (i in 1:(ntt*na)) {
  t = A[i,1]
  a = as.character(A[i,2])
  d = dates[which(dates == dd[[a]][t])]
  anom_t[,,t,a] = d360[,,t,a,"t"] - c360[,,d,"t"]
}

########## Grafico ############################################################
print("Graficando...")
GG = vector("list", ntt*na)
for (i in 1:(ntt*na)) {
  t = A[i,1]
  a = as.character(A[i,2])
  d = dates[which(dates == dd[[a]][t])]
  #print(i);print(t);print(a);print(d)

  df = data.frame(
    lons = rep(lons, ny),
    lats = rep(lats, each = nx),
    #temp = array(d360[,,t,a,"t"] - c360[,,d,"t"]),
    temp = array(anom_t[,,t,a]),
    geop = array(d360[,,t,a,"gh"])
  )
  
  df_hr = data.frame(
    lons_hr = rep(lons_hr, nyhr),
    lats_hr= rep(lats_hr, each = nxhr),
    pv = array(pv_d360[,,t,a])
  )
  
  GG[[i]] = ggplot() +
    geom_contour_fill(
      data = df_hr,
      aes(x = lons_hr, y = lats_hr, z = pv), breaks = brks$pv) +
    scale_fill_gradientn(
      colours = pv_cols,
      breaks = brks$pv, labels = labs$pv, limits = range(brks$pv)
    ) +
    geom_contour2(data = df, aes(x = lons, y = lats, z = geop), breaks = brks$gh, linetype = "solid", size = 0.225, alpha = 0.65) +
    geom_contour2(data = df, aes(x = lons, y = lats, z = temp, color = ..level..), breaks = 1, linetype = "solid", size = 0.25) +
    geom_contour2(data = df, aes(x = lons, y = lats, z = temp, color = ..level..), breaks = brks$t[which(brks$t > 0)], linetype = "solid", size = 0.2) +
    geom_contour2(data = df, aes(x = lons, y = lats, z = temp, color = ..level..), breaks = brks$t[which(brks$t < 0)], linetype = "dashed", size = 0.2) +
    scale_color_gradientn(
      colours = t_cols,
      breaks = brks$t,
      limits = range(brks$t),
      guide = "none"
    ) +
    guides(
      fill = guide_colorstrip(
        title = "q [PVU]     ", barwidth = 14, barheight = 0.3, ticks = TRUE, ticks.color = "gray50",
        label.vjust = 2
      )
    ) +
    geom_map(data = SH, aes(x = long, y = lat, map_id = region), map = SH, fill = NA, color = "black", size = 0.05) +
    coord_map(projection = "stereographic", orientation = c(-90,0,0), ylim = c(-90,lats[1])) +
    scale_y_continuous(breaks = c(-60,-30), name = NULL, labels = NULL) +
    scale_x_continuous(breaks = NULL,name = NULL, labels = NULL) +
    geom_vline(xintercept = seq(-120,180,60), size = 0.05, linetype = "dashed", color = "white") +
    annotate("text", x = -120, y = -30, label = "30°S", color = "white", size = 1.2) +
    annotate("text", x = -120, y = -60, label = "60°S", color = "white", size = 1.2) +
    geom_ribbon(data = df, aes(x = lons, y = lats, ymin = lats[1], ymax = 25), fill = "white") +
    #geom_text_contour(data = df, aes(x = lons, y = lats, z = geop), breaks = 29.5, size = 2.1, stroke = 0.1) +
    geom_text_contour(data = df, aes(x = lons, y = lats, z = temp), breaks = c(20,40,60), size = 2.1, stroke = 0.1, color = "red") +
    #geom_text_contour(data = df, aes(x = lons, y = lats, z = geop), breaks = 30, size = 1.7, stroke = 0.1) +
    labs(title = paste0(ttl[[a]][t], "\n[", stl[t], "]")) +
    theme(
      panel.ontop = TRUE, 
      panel.background = element_blank(), 
      panel.grid = element_line(size = 0.075, linetype = "dashed", color = "white"),
      plot.title = element_text(hjust = 0, vjust = -16, size = 5),
      plot.margin = margin(c(-0.6,-0.18,-0.1,-0.18), unit = "cm"),
      legend.position = "bottom", legend.direction = "horizontal",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5),
      axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()
    )
  print(paste0("Graficado ",a,"-",d))
}

print("Combinando ...")
figure = ggarrange(plotlist = GG, ncol = 4, nrow = 6, common.legend = TRUE, legend = "bottom")
print("Graficado.")
print('Guardando ...')
ggsave(
  here("outputs/ftropos", "horiz_field_10hpa_highres_5.png"),
  figure,  width = 16, height = 22, units = "cm", dpi = 400
)
print('Hecho.')
