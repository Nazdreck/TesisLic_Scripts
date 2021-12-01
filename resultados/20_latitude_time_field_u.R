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

########## Predefiniendo variables ############################################
lats = seq(15,-60, -2.5); ny = length(lats)
plvl = 1
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "2000-01-01")), 6, 10)
stats = c("2002", "2019", "climatology"); ns = length(stats)
years = c(2002,2019); na = length(years)

########## Extrayendo datos ###################################################
DATA = array(double(ny*nt*ns), c(nt,ny,ns), dimnames = list(dates = dates, lats = lats, stats = stats))


for (s in stats) {
  nc = nc_open(
    here("inputs/era5_plvls", paste0("u", "_", s, ".nc"))
  )
  
  data_extracted = ncvar_get(
    nc,
    start = c(
      1,
      which(ncvar_get(nc, "lats") == lats[1]),
      which(ncvar_get(nc, "plvl") == plvl),
      1
    ),
    count = c(
      c(-1,ny,1,-1)
    )
  )
  nc_close(nc); rm(nc)
  
  data_extracted = t(apply(data_extracted, c(2,3), mean))
  
  DATA[,,s] = data_extracted
}

rm(data_extracted)

########## Anomalias ##########################################################
ANOM = array(double(nt*ny*na), c(nt,ny,na), list(dates = dates, lats = lats, years = years))

for (a in as.character(years)) {
  ANOM[,,a] = DATA[,,a] - DATA[,,"climatology"]
}

########## Elementos para graficar ############################################
mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")

brks = seq(-110,110,10); nb = length(brks)
labs = as.character(brks)
labs[which(odd(1:nb))] = " "
brks_p = seq(30,90,30)
brks_n = seq(-60,-30,30)

cols = c(
  sequential_hcl((nb-1)/2, palette = "PuRd", rev = FALSE),
  sequential_hcl((nb-1)/2, palette = "Greens", rev = TRUE)
)


######### Grafico #############################################################
df = data.frame(
  days = rep(days, ny),
  lats = rep(lats, each = nt),
  u2002 = array(DATA[,,"2002"]),
  u2019 = array(DATA[,,"2019"]),
  uclim = array(DATA[,,"climatology"])
)

GG = ggplot(df) +
  aes(x = days, y = lats) +
  geom_contour_fill(aes(z = uclim), breaks = brks, na.fill = FALSE) +
  scale_fill_gradientn(
    colors = cols,
    breaks = brks, labels = labs, limits = range(brks)
  ) +
  geom_contour2(aes(z=uclim), breaks = 0, size = 0.3, color = "black") +
  guides(
    fill = guide_colorstrip(
      title = "m/s   ", barwidth = 11, barheight = 0.3, title.position = "left", title.hjust = 0.5,
      ticks = TRUE, ticks.color = "gray50", label.vjust = 2
    )
  ) +
  geom_contour2(aes(z = u2002), breaks = brks_p, linetype = "solid", size = 0.1, color = "red") +
  geom_contour2(aes(z = u2002), breaks = brks_n, linetype = "dashed", size = 0.1, color = "red") +
  geom_contour2(aes(z = u2002), breaks = 0, linetype = "solid", size = 0.35, color = "red") + 
  geom_contour2(aes(z = u2019), breaks = brks_p, linetype = "solid", size = 0.1, color = "blue") +
  geom_contour2(aes(z = u2019), breaks = brks_n, linetype = "dashed", size = 0.1, color = "blue") +
  geom_contour2(aes(z = u2019), breaks = 0, linetype = "solid", size = 0.35, color = "blue") + 
  scale_x_continuous(
    breaks = mini, labels = mname, limits = c(45,320),
    expand = c(0,0), name = NULL,
  ) +
  scale_y_latitude(
    breaks = seq(-60,15,15)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.1, color = "gray20") +
  labs(title = expression(bar(u)~" (1 hPa)")) +
  theme_light() +
  theme(
    plot.margin = margin(c(0,0,0,0), unit = "cm"),
    plot.title = element_text(size = 7, vjust = -7, hjust = 0.01, margin = margin(0,0,0,0, unit = "cm")),
    legend.margin = margin(c(0.05,0,0,0), unit = "cm"),
    legend.title = element_text(face = "bold", size = 7),
    legend.text = element_text(size = 5),
    legend.position = "bottom",
    legend.box.spacing = unit(0, units = "cm"),
    axis.title = element_blank(),
    axis.text = element_text(size = 5),
  )

ggsave(
  here("outputs/ftropos", paste0("lat_time_field_",plvl,"hpa.png")),
  GG, device = "png", width = 16, height = 6.75, units = "cm", dpi = 200
)











