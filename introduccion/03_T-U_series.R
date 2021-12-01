# Genera las series temporales para T en 60-90°S y U en 60°S para el periodo
# 1979-2019 Tambien grafica la climatologia del periodo 1979-2018 y sin 2002. Se
# plotean los maximos y los minimos de las series, los percentiles 10, 30, 70 y
# 90 se señalan en sombrado. El script funciona de tal manera que genera un
# objeto ggplot para ser llamado desde otros script y combinar con otros
# graficos
########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ggplot2)        # Data Visualisation
require(colorspace)     # Paletas de Colores.
require(abind)          # Combine Multidimensional Arrays
require(ggpubr)         # "ggplot2" Based Publication Ready Plots
require(ncdf4)          # For manipulate .nc files


########## Predefiniendo variables ############################################
vars = c("t","u"); nv = length(vars)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "1900-01-01")), 6,10)
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = 1979:2019; na = length(years)

DATA = array(double(nt*na*nv), dim = c(nt,na,nv), dimnames = list(days = days, years = years, vars = vars))

########## Extrayendo datos ###################################################
print("Extrayendo datos de T y U...")
DIR = "../../../datos/ERA5_aux/nuevos/daily"
partial = NULL

for (v in vars) {
  for (a in years) {
    for (m in months) {
      if (m == "02") s = 28 else s = -1
      if (v == "t") {lt = -1; pl = 30} else {lt = 1; pl = 10}
      nc = nc_open(
        here(DIR, paste(v, m, a, "day.nc", sep = "_"))
      )
      data_extracted = ncvar_get(
        nc,
        varid = names(nc$var)[2],
        start = c(
          1,
          which(ncvar_get(nc, "latitude") == -60),
          which(ncvar_get(nc, "level") == pl),
          1
        ),
        count = c(-1,lt,1,s)
      )
      nc_close(nc); rm(nc)
      
      if (v == "t") partial = abind(partial, apply(data_extracted, c(3), mean), along = 1)
      if (v == "u") partial = abind(partial, apply(data_extracted, c(2), mean), along = 1)
    }
    DATA[,as.character(a),v] = partial
    partial = NULL
  }
}

rm(data_extracted,partial)
print("Extraido")
########## Climatologia y percentiles #########################################
print("Calculando estadisticos...")
DATA[,,"t"] = DATA[,,"t"]-273
CLIM = apply(DATA[,as.character(c(1979:2001, 2003:2018)),], c(1,3), mean)
MAX = apply(DATA[,as.character(c(1979:2001, 2003:2018)),], c(1,3), max)
MIN = apply(DATA[,as.character(c(1979:2001, 2003:2018)),], c(1,3), min)

pers = c(10,30,70,90); np = length(pers)

PER = array(double(np*nt*nv), c(nt,np,nv), dimnames = list(days = days,percentile = pers, vars = vars))

for (v in vars) {
  for (d in days) {
    for (p in pers) {
      x1 = which(Percentile(DATA[d,,v]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
      x2 = min(DATA[d,,v][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
      PER[d,as.character(p),v] = x2 # Guarda el valor minimo encontrado
    }
  }
}

print("Calculados.")
########## Figuras ############################################################
print("Graficando series temporales...")
lims = list(t = c(-90,-20), u = c(-25,105))
brks = list(t = seq(-90,-20,10), u = seq(-20,100,20))
leg = list(t = "none", u = c(0.015, 0.86))
tit = list(t = expression(bold(bar(T)~ "[°C] (60-90°S, 30 hPa)")), u = expression(bold(bar(u)~ "[m/s] (60°S, 10 hPa)")))
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")
axsx = list(t = element_text(), u = element_blank())


gg_series = vector("list",2); names(gg_series) = vars

for (v in vars) {
  df_lines = data.frame(
    days = rep(days,9),
    vars = c(
      p10 = PER[,"10",v],
      p90 = PER[,"90",v],
      p30 = PER[,"30",v],
      p70 = PER[,"70",v],
      max = MAX[,v],
      min = MIN[,v],
      clim = CLIM[,v],
      d2002 = DATA[,"2002",v],
      d2019 = DATA[,"2019",v]
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
    pmin = c(p10 = PER[,"10",v], p30 = PER[,"30",v]),
    pmax = c(p90 = PER[,"90",v], p70 = PER[,"70",v]),
    group = rep(LETTERS[1:2], each = nt)
  )
  
  gg_series[[v]] = ggplot() +
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
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.1) +
    scale_x_continuous(name = NULL, limits = range(1:365), breaks = mini, labels = mname, expand = c(0,0))+
    scale_y_continuous(
      name = NULL, breaks = brks[[v]], limits = lims[[v]],expand = c(0,0)
    ) +
    labs(
      title = tit[[v]]
    ) +
    theme_light() +
    theme(
      axis.text.x = axsx[[v]],
      legend.key.size = unit(0.275, "cm"),
      legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
      legend.margin = margin(-2, 3, 3, 3),
      legend.key = element_rect(fill = alpha("orange", 0.05)),
      legend.title = element_blank(),
      legend.position = leg[[v]],
      legend.justification = c("left", "top"),
      plot.title = element_text(face = "bold", vjust = -7, hjust = 0.02, size = 7, margin = margin(t = -10)),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6),
    )
}

print("Graficado.")

fig_series = ggarrange(
  plotlist = list(gg_series[[2]],gg_series[[1]]), ncol = 1, nrow = 2, align = "v"
)

ggsave(
  here("outputs/introduccion", "T-U_series.png"),
  fig_series, device = "png", width = 12, height = 9, units = "cm", dpi = 200
)

