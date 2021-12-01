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

########## Extraccion de datos ################################################
DIR = "../../../datos/ERA5_aux/nuevos/daily"

lons = seq(0,357.5,2.5); nx = length(lons)
years = 1979:2019; na = length(years)
months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
days = 1:365; nt = length(days)
dates = as.character(as.Date(days, origin = "2018-12-31"))

DATA = array(double(nx*nt*na), dim = c(nx,nt,na), 
             dimnames = list(lons = lons, dates = dates, years = years)
)
partial = NULL

for (a in years) {
  for (m in months) {
    if (m == "02") s = 28 else s = -1
    nc = nc_open(
      here(DIR, paste("gh",m,a, "day.nc", sep = "_"))
    )
    data_extracted = ncvar_get(
      nc,
      varid = names(nc$var)[2],
      start = c(
        1,
        which(ncvar_get(nc, "latitude") == -60),
        which(ncvar_get(nc, "level") == 10),
        1
      ),
      count = c(
        -1,1,1,s
      )
    )
    nc_close(nc); rm(nc)
    
    partial = abind(partial, data_extracted, along = 2)
  }
  DATA[,,as.character(a)] = partial
  partial = NULL
}

# Carga los dias con minimos de flujo de calor
load("inputs/heatflux_mindays.RData")

########## Analisis Espectral ###############################################
GH = DATA/9800
dimnames(GH) = list(lons = lons, dates = dates, years = years)

spectral = array(double(2*nt*2*na), dim = c(nt,2,2,na), 
                 dimnames = list(dates = dates, c("amp", "phase"), 
                                 wavenumber = c(1,2), years = years)
)


for (k in 1:2) {
  for (j in 1:na) {
    for (i in 1:nt) {
      fourier = FitWave(GH[,i,j], k = k)
      spectral[i, , k, j] = c(fourier$amplitude, fourier$phase)
    }
  }
}


########## Percentiles ########################################################
pers = c(10,30,70,90); np = length(pers)

PER = array(double(np*nt*2), c(nt,np,2), dimnames = list(days = days,percentile = pers, k = c(1,2)))

for (k in 1:2) {
  for (d in 1:nt) {
    for (p in pers) {
      x1 = which(Percentile(spectral[d,"amp", k,]) > p/100) # Ubicaciones de los elementos por encima del percentil p%
      x2 = min(spectral[d,"amp", k,][x1]) # Busca el valor mas chico de los elementos por encima del percentil p%
      PER[d,as.character(p),k] = x2 # Guarda el valor minimo encontrado
    }
  }
}


########## Figuras ############################################################
mini = c(1,32,63,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")

########## Graficos de amplitud

brks = list(seq(0,2.25,0.25), seq(0,1.25,0.25))
lims = list(c(-0.125,2.25), c(-0.0625,1.25))
leg = list(c(0.015, 0.88), "none")
tit = list(expression(Z[1]~ "[km] (60°S, 10 hPa)"), expression(Z[2]~ "[km] (60°S, 10 hPa)"))
axsx = list(element_blank(), element_text())

GG_amp = vector("list", 2)

for (k in 1:2) {
  
  df_lines = data.frame(
    days = rep(days, 9),
    vars = c(
      p10 = PER[,"10",k],
      p90 = PER[,"90",k],
      p30 = PER[,"30",k],
      p70 = PER[,"70",k],
      max = apply(spectral[,"amp",k,], 1, max),
      min = apply(spectral[,"amp",k,], 1, min),
      clim = apply(spectral[,"amp",k,], 1, mean),
      d2002 = spectral[,"amp",k,"2002"],
      d2019 = spectral[,"amp",k,"2019"]
    ),
    group = rep(LETTERS[1:9], each = nt),
    cols = c(
      rep("10-90%", each = nt*2), rep("30-70%", each = nt*2), rep("max/min", each = nt*2),
      rep("clim", each = nt), rep("2002", each = nt), rep("2019", each = nt)
    )
  )
  
  df_lines$cols = factor(df_lines$cols, levels = c("10-90%", "30-70%", "max/min", "clim", "2002", "2019"))
  
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
      name = "Leyenda",
      values = c(
        "max/min" = "darkgreen", "clim" = "black", 
        "2019" = "blue", "2002"= "red", 
        "10-90%"= "gray70", "30-70%" = "gray30"
      )
    ) +
    scale_size_manual(
      values = c(
        "max/min" = 0.15, "clim" = 0.3, 
        "2019" = 0.3, "2002"= 0.3, 
        "10-90%"= 0.1, "30-70%" = 0.1
      )
    ) +
    guides(size = "none") +
    scale_x_continuous(name = NULL, limits = range(77:350), breaks = mini, labels = mname, expand = c(0,0))+
    scale_y_continuous(
      name = NULL, breaks = brks[[k]], limits = lims[[k]],
      expand = c(0,0)
    ) +
    geom_vline(xintercept = mins2002, color = "red", linetype = "dashed", size = 0.1) +
    geom_vline(xintercept = mins2019, color = "blue", linetype = "dashed", size = 0.1) +
    labs(
      title = tit[[k]]
    ) +
    theme_light() +
    theme(
      axis.text.x = axsx[[k]],
      legend.key.size = unit(0.275, "cm"),
      legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
      legend.margin = margin(-2, 3, 3, 3),
      legend.key = element_rect(fill = alpha("orange", 0.05)),
      legend.title = element_blank(),
      legend.position = leg[[k]],
      legend.justification = c("left", "top"),
      plot.title = element_text(face = "bold", vjust = -7, hjust = 0.02, size = 7, margin = margin(t = -10)),
      axis.text = element_text(size = 6),
      legend.text = element_text(size = 6),
    )
}

fig_amp = ggarrange(
  plotlist = list(GG_amp[[1]],GG_amp[[2]]), ncol = 1, nrow = 2
)

ggsave(
  here("outputs/ftropos", "spectral_amplitude2.png"),
  fig_amp, device = "png", width = 12, height = 9, units = "cm", dpi = 200
)

########## Grafico de fase

df_phase = data.frame(
  days = rep(days, 4),
  vars = c(
    spectral[,"phase", 1, "2002"],
    spectral[,"phase", 2, "2002"],
    spectral[,"phase", 1, "2019"],
    spectral[,"phase", 2, "2019"]
  )*180/pi,
  group = rep(LETTERS[1:4], each = nt),
  cols = rep(c("2002","2019"), each = 2*nt),
  ltype = rep(c("k1", "k2", "k1", "k2"), each = nt)
)

GG_phase = ggplot(df_phase) +
  aes(x = days) +
  geom_line(aes(y = vars, group = group, color = cols, linetype = ltype), size = 0.2) +
  facet_grid(rows = vars(cols), scales = "free") +
  scale_color_manual(values = c("2002" = "red", "2019" = "blue")) +
  scale_linetype_discrete(name = "N° Onda\nPlanetario:", labels = c("k = 1", "k = 2")) +
  guides(color = "none") +
  labs(title = "\u03C6 [°E] (60°S, 10 hPa)") +
  scale_x_continuous(name = NULL, limits = range(136:289), breaks = mini, labels = mname, expand = c(0,0)) +
  scale_y_continuous(name = NULL, limits = c(0,360), breaks = seq(0,360,60), labels = paste0(seq(0,360,60),"°"), expand = c(0,0)) +
  theme_light() +
  theme(
    legend.key.size = unit(0.35, "cm"),
    legend.background = element_rect(fill = alpha("orange", 0.5), color = "black", size = 0.1),
    legend.margin = margin(3, 3, 3, 3),
    legend.key = element_rect(fill = alpha("orange", 0.01)),
    legend.position = c(0.99,0.99),
    legend.direction = "horizontal",
    legend.justification = c("right", "top"),
    plot.title = element_text(face = "bold", margin = margin(t = 10, b = -10), hjust = 0.02, size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    strip.text = element_text(face = "bold", color = "black", size = 6, margin = margin(r = 0, l = 0))
  )

ggsave(
  here("outputs/ftropos", "spectral_phase.png"),
  GG_phase, device = "png", width = 12, height = 6.75, units = "cm", dpi = 200
)
