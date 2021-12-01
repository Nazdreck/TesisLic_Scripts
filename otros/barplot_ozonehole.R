# Genera un bar plot para el area promedio del agujero de ozono entre septiembre
# y octubre, para cada año del periodo 1979-2019. EL script funciona de tal
# manera que genera un objeto ggplot para ser llamado desde otros script y
# combinar con otros graficos
########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ggplot2)        # Data Visualisation
require(colorspace)     # Paletas de Colores.
require(ncdf4)          # For manipulate .nc files


source(here("functions/area_on_SH.R"))

########## Predefiniendo las variables #######################################
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(0,-90,-2.5); ny = length(lats)
days = 1:61; nd = length(days)
years = 1979:2019; na = length(years)

########## Extrayendo datos ##################################################
print("Extrayendo datos de ozono...")
nc_ozone = nc_open(
  here("inputs/era5_ozone/tco_09-10_79-19.nc")
)

DATA = ncvar_get(
  nc_ozone,
  varid = "tco",
  start = c(
    1,
    which(ncvar_get(nc_ozone, "lats") == 0),
    1,1
  ),
  count = c(-1,-1,-1,-1)
)

nc_close(nc_ozone)
rm(nc_ozone)

dimnames(DATA) = list(lons = lons, lats = lats, days = days, years = years)
print("Extraido.")
########## Calculo del area del agujero de ozono ##############################
print("Calculando area de agujero de ozono...")
OZONE = array(double(nx*ny*nd*na), dim = c(nx,ny,nd,na), dimnames = list(lons = lons, lats = lats, days = days, years = years))

hole = which(DATA/2.1415e-5 < 220, arr.ind = TRUE)

n = nrow(hole)
for (i in 1:n){
  cd = hole[i,]
  OZONE[cd[1],cd[2],cd[3],cd[4]] = 1
}

area = array(double(nd*na), dim = c(nd,na), dimnames = list(days=days, years=years))

for (a in 1:na){
  for (d in 1:nd){
    area[d,a] = area.on.SH(OZONE[,,d,a],2.5)
  }
}

area_mean = apply(area, 2, mean)

rm(list = c("OZONE","DATA","hole"))
print("Calculado.")
########## Figura #############################################################
print("Creando un grafico de barras...")
df = data.frame(
  years = years,
  area = area_mean/1e6,
  mark = c(rep("A", 23),"B", rep("A",16), "C"),
  rank = rank(area_mean) - 1
)

gg_ozone = ggplot(df) +
  aes(
    x = years,
    y = area,
    fill = mark
  ) +
  geom_col(width = 0.7, show.legend = FALSE, color = "black", size = 0.15) +
  scale_fill_manual(values = alpha(c("gray","red","blue"), 0.8)) +
  geom_text(aes(label = rank), vjust = -0.25, size = 1.7) +
  geom_hline(yintercept = 17.84, linetype = "dashed", size = 0.2) +
  annotate("text", y = 17.84, x = 1981, label = "Superficie\nSudamérica", size = 1.8) +
  scale_x_continuous(
    breaks = c(1979,seq(1983,2015,4),2019),
    expand = c(0,0),
    limits = c(1978,2020),
    name = NULL
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = seq(0,25,5),
    limits = c(0,28),
    expand = c(0,0)
  ) +
  labs(title = expression("Área OZ [millones de km"^{2}*"]")) +
  theme_light() +
  theme(
    plot.title = element_text(face = "bold", vjust = -7, hjust = 0.02, size = 7, margin = margin(t = -10)),
    axis.text = element_text(size = 6)
  )
print("Creado.")


ggsave(
  here("outputs/introduccion", "TCO_barplot.png"),
  gg_ozone, device = "png", width = 8, height = 4.5, units = "cm", dpi = 200
)


