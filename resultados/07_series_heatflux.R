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

source(here("functions/get_armonic.R"))

########## Predefiniendo variables ############################################
vars = c("t", "v"); nv = length(vars)
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(-45,-75,-2.5); ny = length(lats)
days = 1:365; nt = length(days)
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = 1979:2019; na = length(years)

DIR = "../../../datos/ERA5_aux"

DATA = array(double(nx*ny*nt*na*nv), dim = c(nx,ny,nt,na,nv), dimnames = list(lons = lons, lats = lats, days = days, years = years, vars = vars))

########## Extraccion de datos ################################################
partial = NULL

for (v in vars) {
  for (a in years) {
    for (m in months) {
      if (m == "02") s = 28 else s = -1
      nc = nc_open(
        here(DIR, paste(v, m, a, "day.nc", sep = "_"))
      )
      data_extracted = ncvar_get(
        nc,
        varid = names(nc$var)[2],
        start = c(
          1,
          which(ncvar_get(nc, "latitude") == -45),
          which(ncvar_get(nc, "level") == 10),
          1
        ),
        count = c(
          -1, ny, 1, s
        )
      )
      nc_close(nc); rm(nc)
      
      partial = abind(partial, data_extracted, along = 3)
    }
    DATA[,,,as.character(a),v] = partial
    partial = NULL
  }
}

rm(partial); rm(data_extracted)

########## Flujo de calor #####################################################
data_t = DATA[,,,,"t"]
data_v = DATA[,,,,"v"]

hf = array(double(ny*nt*na), dim = c(ny,nt,na), dimnames = list(lats = lats, days = days, years = years))

for (a in 1:na) {
  for (t in 1:nt) {
    for (y in 1:ny) {
      hf[y,t,a] = cov(data_t[,y,t,a],data_v[,y,t,a])
    }
  }
}

hf = apply(hf, c(2,3), mean)

########## Percentiles ########################################################
pers = c(10,30,70,90); np = length(pers)

PER = array(double(np*nt), c(nt,np), dimnames = list(days = days,percentile = pers))


for (d in days) {
  for (p in pers) {
    x1 = which(Percentile(hf[d,]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
    x2 = min(hf[d,][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
    PER[d,as.character(p)] = x2 # Guarda el valor minimo encontrado
  }
}

########## Minimos locales ####################################################
mins2002 = as.numeric(names(which(hf[which(diff(sign(diff(hf[,"2002"]))) == 2)+1,"2002"] < -100))); nm2 = length(mins2002); names(mins2002) = NULL
mins2019 = as.numeric(names(which(hf[which(diff(sign(diff(hf[,"2019"]))) == 2)+1,"2019"] < -100))); nm9 = length(mins2019); names(mins2019) = NULL

save(list = c("mins2002","mins2019"),file = "inputs/heatflux_mindays.RData")
######## Grafico ##############################################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
brks = seq(-400,100,100)


df_lines = data.frame(
  days = rep(days, 9),
  vars = c(
    p10 = PER[,"10"],
    p90 = PER[,"90"],
    p30 = PER[,"30"],
    p70 = PER[,"70"],
    max = apply(hf[,-c(24,41)], 1, max),
    min = apply(hf[,-c(24,41)], 1, min),
    clim = apply(hf[,-c(24,41)], 1, mean),
    d2002 = hf[,"2002"],
    d2019 = hf[,"2019"]
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

df_points = data.frame(
  days = c(mins2002,mins2019),
  hf = c(hf[mins2002,"2002"],hf[mins2019,"2019"]),
  year = c(rep("2002", nm2), rep("2019", nm9))
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
  geom_point(
    data = df_points, aes(x = days, y = hf, group = year, color = year), show.legend = FALSE, shape = 4, size = 1
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
  scale_x_continuous(name = NULL, limits = range(91:335), breaks = mini, labels = mname, expand = c(0,0))+
  scale_y_continuous(
    name = NULL, breaks = brks, limits = range(brks),
    expand = c(0,0)
  ) +
  labs(
    title = expression(bar("v'T'")~ "[K ms"^{-1}*"] (45-75°S, 10 hPa)")
  ) +
  theme_light() +
  theme(
    legend.key.size = unit(0.275, "cm"),
    legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
    legend.margin = margin(-2, 3, 3, 3),
    legend.key = element_rect(fill = alpha("orange", 0.05)),
    legend.title = element_blank(),
    legend.position = c(0.01,0.025),
    legend.justification = c("left", "bottom"),
    plot.title = element_text(face = "bold", vjust = -7.5, hjust = 0.02, size = 7, margin = margin(t = -10)),
    axis.text = element_text(size = 6),
    legend.text = element_text(size = 6),
  )


ggsave(
  here("outputs/ftropos", "hf_series.png"),
  GG, device = "png", width = 12, height = 5.0625, units = "cm", dpi = 200
)










