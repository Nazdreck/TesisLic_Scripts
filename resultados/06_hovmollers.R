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
lons = seq(0,357.5,2.5); nx = length(lons)
plvl = c(500,10); nz = length(plvl)
days = 1:123; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "1900-07-01")), 6,10)
dfull = substr(as.character(as.Date(1:365-1, origin = "1900-01-01")), 6,10)
years = c("2002", "2019", "climatology"); na = length(years)

DATA = array(
  double(nx*nz*na*nt), c(nx,nz,nt,na), 
  dimnames = list(lons = lons, plvl = plvl, dates = dates, years = years)
)

########## Extrayendo datos ###################################################

for (a in years) {
  nc = nc_open(
    here("inputs/era5_plvls", paste0("gh_",a,".nc"))
  )
  for (z in plvl) {
    data_extracted = ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == -60),
        which(ncvar_get(nc, "plvl") == z),
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "07-01")
      ),
      count = c(
        -1,1,1,nt
      )
    )
    
    DATA[,as.character(z),,a] = data_extracted
  }
  nc_close(nc); rm(nc)
}
rm(data_extracted)


# Carga los dias con minimos de flujo de calor
load("inputs/heatflux_mindays.RData")

########## Analisis Espectral #################################################
waves = c("k1", "k2"); nw = length(waves)

FOURIER = array(
  double(nx*nz*nt*na*nw), c(nx,nz,nt,na,nw),
  dimnames = list(lons = lons, plvl = plvl, days = days, years = years, waves = waves)
)

FOURIER[,,,,"k1"] = apply(DATA/9.8, 2:4, FilterWave, k = 1)
FOURIER[,,,,"k2"] = apply(DATA/9.8, 2:4, FilterWave, k = 2)
########## Elementos para graficar
dnumb = c(1,15,32,46,63,77,93,107,123)
dname = c("1-JUL","15-JUL","1-AGO","15-AGO", "1-SEP", "15-SEP", "1-OCT", "15-OCT", "31-OCT")
hfmins = list("2002" = (mins2002-181)[which(mins2002-181 > 0)], "2019" = (mins2019-181)[which(mins2019-181 > 0)], "climatology" = c())

# BREAKS
brks = list(
  "10" = seq(-2500,2500,500),
  "500" = seq(-300,300,50)
)
nb = list("10" = length(brks$`10`), "500" = length(brks$`500`))

# LABELS
labs = list(
  "10" = c(-1000,-500,0,500,1000),
  "500" = c(-200,-100,0,100,200)
)


# COLORES
cols = list(
  "10" = c(
    sequential_hcl((nb$`10`-1)/2, palette = "Blues", rev = TRUE),
    sequential_hcl((nb$`10`-1)/2, palette = "Reds")
  ),
  "500" = c(
    sequential_hcl((nb$`500`-1)/2, palette = "Greens", rev = TRUE),
    sequential_hcl((nb$`500`-1)/2, palette = "Oranges")
  )
)

# TITULO
ttl = c("2002", "2019", "Climatología")

# Matriz iteradora
# A = matrix(
#   c(1, 10, 1, 1, 500, 1, 2, 10, 1, 2, 500, 1,
#     1, 10, 2, 1, 500, 2, 2, 10, 2, 2, 500, 2,
#     1, 10, 3, 1, 500, 3, 2, 10, 3, 2, 500, 3),
#   ncol = 3, nrow = 12, byrow = TRUE
# )

A = matrix(
  c(1, 10, 1, 2, 10, 1, 1, 500, 1, 2, 500, 1,
    1, 10, 2, 2, 10, 2, 1, 500, 2, 2, 500, 2,
    1, 10, 3, 2, 10, 3, 1, 500, 3, 2, 500, 3),
  ncol = 3, nrow = 12, byrow = TRUE
)

########## Graficos
GG = vector("list", nz*nw*na)
for (i in 1:(nz*nw*na)) {
  a = A[i,3]
  z = as.character(A[i,2])
  k = A[i,1]
  
  df = data.frame(
    lons = rep(ConvertLongitude(c(lons[73:144],lons[1:72])), nt),
    days = rep(days, each = nx),
    geop = array(abind(FOURIER[73:144,z,,a,k],FOURIER[1:72,z,,a,k], along = 1))
  )
  
  GG[[i]] = ggplot(df) +
    aes(x = lons, y = days, z = geop) +
    geom_contour_fill(breaks = brks[[z]]) +
    geom_text_contour(breaks = labs[[z]], stroke = 0.05, size = 1.8) +
    scale_fill_gradientn(
      colors = cols[[z]],
      breaks = brks[[z]],
      limits = range(brks[[z]])
    ) +
    guides(
      fill = "none"
    ) +
    scale_x_longitude(breaks = c(-120,-60,0,60,120)) +
    scale_y_continuous(
      name = NULL, breaks = dnumb, labels = dname, expand = c(0,0),
      trans = "reverse"
    ) +
    geom_hline(yintercept = hfmins[[a]], color = "white", linetype = "dashed", size = 0.2) +
    labs(title = bquote(Z[.(k)]~ " ("~.(ttl[a])~") "~ .(z)~ "hPa")) +
    theme_light() +
    theme(
      # panel.ontop = TRUE, 
      # panel.background = element_blank(), 
      # panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray85"),
      # panel.grid.minor = element_blank(),
      plot.margin = margin(c(0,0,-0.45,0), unit = "cm"),
      plot.title = element_text(hjust = 0, vjust = -5, size = 5),
      axis.text = element_text(size = 4)
    )
  
}

figure = ggarrange(plotlist = GG, ncol = 4, nrow = 3)

print('Guardando ...')
ggsave(
  here("outputs/ftropos", "hovmollers2.png"),
  figure,  width = 16, height = 22, units = "cm", dpi = 200
)
print('Hecho.')

