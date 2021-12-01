########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ggplot2)        # Data Visualisation
require(colorspace)     # Paletas de Colores.
require(ncdf4)          # For manipulate .nc files


source(here("functions/area_on_SH.R"))

########## Predefiniendo variables ############################################
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(0,-90,-2.5); ny = length(lats)
plvl = 10  
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-01-01")),6,10)
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = 1979:2019; na = length(years)

DATA = array(
  double(nx*ny*nt*na), dim = c(nx,ny,nt,na), 
  dimnames = list(lons = lons, lats = lats, dates = dates, years = years)
)
########## Extraccion de datos ################################################
DIR = "../../../datos/ERA5_aux/nuevos/daily"
partial = NULL

for (a in years) {
  for (m in months) {
    if (m == "02") s = 28 else s = -1
    nc = nc_open(
      here(DIR, paste("pv", m, a, "day.nc", sep = "_"))
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
      count = c(-1,ny,1,s)
    )
    nc_close(nc); rm(nc)
    
    partial = abind(partial, data_extracted, along = 3)
  }
  DATA[,,,as.character(a)] = partial
  partial = NULL
}


rm(data_extracted,partial)
print("Extraido")
########## Calculo de area del vortice ########################################
VORTEX = array(
  double(nx*ny*nt*na), dim = c(nx,ny,nt,na), 
  dimnames = list(lons = lons, lats = lats, dates = dates, years = years)
)

# Ubica los puntos donde estan los valores menores a -400 PVU (en tiempo y espacio)
v400 = which(DATA/10^-6 < -400, arr.ind = TRUE)

n = nrow(v400) # Cantidad total de puntos
for (i in 1:n){
  cd = v400[i,] # Toma la ubicacion de cada punto (en tiempo y espacio)
  VORTEX[cd[1],cd[2],cd[3],cd[4]] = 1 # Indexa y asigna 1 donde se cumple el criterio
}

AREA = array(double(nt*na), dim = c(nt,na), dimnames = list(dates=dates, years=years))

for (a in 1:na){
  for (t in 1:nt){
    AREA[t,a] = area.on.SH(VORTEX[,,t,a],2.5)
  }
}
########## Promedio movil #####################################################
sds = 2 # Promedio 5 dias, centrado.

for (a in 1:na) {
  for (i in (sds+1):(nt-sds)) {
    AREA[i, a] = mean(AREA[(i-sds):(i+sds), a])
  }
}

AREA[c(1:sds,(nt-sds+1):nt),] = NA
AREA = AREA/1e6
########## Percentiles ########################################################
pers = c(10,30,70,90); np = length(pers)
PER = array(double(np*nt), c(nt,np), dimnames = list(days = days,percentile = pers))

for (t in days[(sds+1):(nt-sds)]) {
  for (p in pers) {
    x1 = which(Percentile(AREA[t,]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
    x2 = min(AREA[t,][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
    PER[t,as.character(p)] = x2 # Guarda el valor minimo encontrado
  }
}

PER[c(1:sds,(nt-sds+1):nt),] = NA

########## Figuras ############################################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
aa = as.character(c(1979:2001, 2003:2018))

brks = seq(0,60,10)

df_lines = data.frame(
  days = rep(days, 9),
  vars = c(
    p10 = PER[,"10"],
    p90 = PER[,"90"],
    p30 = PER[,"30"],
    p70 = PER[,"70"],
    max = apply(AREA[,aa], 1, max),
    min = apply(AREA[,aa], 1, min),
    clim = apply(AREA[,aa], 1, mean),
    d2002 = AREA[,"2002"],
    d2019 = AREA[,"2019"]
  ),
  group = rep(LETTERS[1:9], each = nt),
  cols = c(
    rep("10-90%", each = nt*2), rep("30-70%", each = nt*2), rep("Máx/mín", each = nt*2),
    rep("Clim.", each = nt), rep("2002", each = nt), rep("2019", each = nt)
  )
)

df_lines$cols = factor(df_lines$cols, levels = c("10-90%", "30-70%", "Máx/mín", "Clim.", "2002", "2019"))

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
    name = "Leyenda",
    values = c(
      "Máx/mín" = "darkgreen", "Clim." = "black", 
      "2019" = "blue", "2002"= "red", 
      "10-90%"= "gray70", "30-70%" = "gray30"
    )
  ) +
  scale_size_manual(
    values = c(
      "Máx/mín" = 0.15, "Clim." = 0.3, 
      "2019" = 0.3, "2002"= 0.3, 
      "10-90%"= 0.1, "30-70%" = 0.1
    )
  ) +
  guides(size = "none") +
  scale_x_continuous(name = NULL, limits = range(120:336), breaks = mini, labels = mname, expand = c(0,0))+
  scale_y_continuous(
    name = NULL, breaks = brks, limits = range(brks),
    expand = c(0,0)
  ) +
  geom_hline(yintercept = 17.84, linetype = "dashed", size = 0.2) +
  labs(
    title = expression("Área VEP [millones de km"^{2}*"]")
  ) +
  theme_light() +
  theme(
    legend.key.size = unit(0.25, "cm"),
    legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
    legend.margin = margin(-2, 3, 3, 3),
    legend.key = element_rect(fill = alpha("orange", 0.05)),
    legend.title = element_blank(),
    legend.position = c(0.01,0.90),
    legend.justification = c("left", "top"),
    plot.title = element_text(face = "bold", vjust = -7.5, hjust = 0.02, size = 7, margin = margin(t = -10)),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 5),
  )


ggsave(
  here("outputs/ftropos", "vortex_area.png"),
  GG, device = "png", width = 12, height = 5.0625, units = "cm", dpi = 200
)






















