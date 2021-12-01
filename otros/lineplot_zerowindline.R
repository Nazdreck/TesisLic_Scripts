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

source(here("functions","is_leapyear.R"))   # Ask if year is leap

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(10,-90,-2.5); ny = length(lats)
plvl = c(100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
days = 1:92; nt = length(days)
dates = substr(as.character(as.Date(days - 1, origin = "1900-07-01")),6,10)
months = c("07","08","09") 
years = 1979:2019; na = length(years)
quin = c("07-1","07-2","08-1","08-2","09-1","09-2"); nq = length(quin)

DATA = array(
  double(ny*nz*nq*na), c(ny,nz,nq,na),
  dimnames = list(lats = lats, plvl = plvl, quin = quin, years = years)
)

########## Extrayendo datos ###################################################
print("Extrayendo datos de T y U...")
DIR = "../../../datos/ERA5_aux/nuevos/daily"
partial = NULL


for (a in years) {
  for (m in months) {
    nc = nc_open(
      here(DIR, paste("u", m, a, "day.nc", sep = "_"))
    )
    data_extracted = ncvar_get(
      nc,
      varid = names(nc$var)[2],
      start = c(
        1,
        which(ncvar_get(nc, "latitude") == lats[1]),
        which(ncvar_get(nc, "level") == plvl[1]),
        1
      ),
      count = c(-1,-1,-1,-1)
    )
    nc_close(nc); rm(nc)
    
    prom = abind(
      apply(data_extracted[,,,1:15], c(2,3), mean),
      apply(data_extracted[,,,-(1:15)], c(2,3), mean),
      along = 3
    )
    
    partial = abind(partial, prom, along = 3)
  }
  DATA[,,,as.character(a)] = partial
  partial = NULL
}


rm(data_extracted,partial,prom)
print("Extraido")

########## QBO ################################################################
QBO = array(vector("character", nq*na), c(nq,na), dimnames = list(quin = quin, years = years))

for (a in as.character(years)) {
  QBO30_E = which(DATA["0","30",,a] > 0, arr.ind = TRUE)
  QBO30_W = which(DATA["0","30",,a] < 0, arr.ind = TRUE)
  QBO[QBO30_E,a] = "QBO30-E"
  QBO[QBO30_W,a] = "QBO30-W"
}

QBO[,"2002"] = "2002"
QBO[,"2019"] = "2019"

########## Elementos para graficos ############################################
tit = c("1-15/JUL", "16-31/JUL", "1-15/AGO", "16-31/AGO", "1-15/SEP", "16-30/SEP")

##### theming
# x axis
tm.axis_text_x_bottom = list("top" = element_text(), "bot" = element_blank())
tm.axis_ticks_x_bottom = list("top" = element_line(), "bot" = element_blank())
tm.axis_ticks_x_top = list("top" = element_blank(), "bot" = element_line())
# y axis
tm.axis_text_y_left = list("left" = element_text(), "center" = element_blank(), "right" = element_blank())
tm.axis_text_y_right = list("left" = element_blank(), "center" = element_blank(), "right" = element_text())
tm.axis_ticks_y_left = list("left" = element_line(), "center" = element_blank(), "right" = element_blank())
tm.axis_ticks_y_right = list("left" = element_blank(), "center" = element_blank(), "right" = element_line())
########## Figura #############################################################

A = matrix(
  c("left", "top", "center", "top", "right", "top", "left", "bot", "center", "bot", "right", "bot"),
  ncol = 2, nrow = 6, byrow = TRUE
)

GG = vector("list", nq); names(GG) = quin

for (q in 1:nq) {
  v = A[q,2]
  h = A[q,1]
  
  df = data.frame(
    lats = rep(rep(lats,nz), na),
    plvl = rep(rep(plvl, each = ny), na),
    uwind = array(DATA[,,q,]),
    group = rep(as.character(years), each = nz*ny),
    colors = rep(QBO[q,], each = nz*ny)
  )
  
  GG[[q]] = ggplot() +
    geom_contour2(
      data = df, aes(x = lats, y = plvl, z = uwind, group = group, color = colors, size = colors),
      breaks = 0
    ) +
    scale_color_manual(
      name = "Leyenda",
      values = c(
        "2019" = "blue", "2002"= "red",
        "QBO30-E" = "orange", "QBO30-W" = "purple"
      )
    ) +
    scale_size_manual(
      values = c(
        "2019" = 0.2, "2002"= 0.2,
        "QBO30-E" = 0.1, "QBO30-W" = 0.1
      )
    ) +
    guides(
      size = "none"
    ) +
    scale_y_level(
      name = NULL, breaks = c(100,70,50,30,20,10,7,5,3,2,1), sec.axis = dup_axis(name = NULL)
    ) +
    scale_x_latitude(
      breaks = seq(-60,0,20), limits = c(-70,10), sec.axis = dup_axis(name = NULL)
    ) +
    annotate("text", x = -60, y = 82, label = tit[q], color = "black", size = 2) +
    annotate("point", x = c(-60,0), y = c(10,30), size = 1) +
    theme_light() +
    theme(
      panel.ontop = TRUE, 
      panel.background = element_blank(), 
      panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray50"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
      plot.title = element_blank(),
      legend.title = element_text(size = 5, face = "bold"),
      legend.text = element_text(size = 4),
      legend.position = "bottom",
      legend.margin = margin(c(0,0,0,0), unit = "cm"),
      legend.box.spacing = unit(0, units = "cm"),
      axis.title = element_blank(),
      # x-axis
      axis.text.x.bottom = tm.axis_text_x_bottom[[v]],
      axis.text.x.top = element_blank(),
      axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[v]],
      axis.ticks.x.top = tm.axis_ticks_x_top[[v]],
      # y-axis
      axis.text.y.left = tm.axis_text_y_left[[h]],
      axis.text.y.right = tm.axis_text_y_right[[h]],
      axis.ticks.y.left = tm.axis_ticks_y_left[[h]],
      axis.ticks.y.right = tm.axis_ticks_y_right[[h]],
      axis.text = element_text(size = 5)
    )
}

figure = ggarrange(plotlist = GG, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/ftropos/zero_wind_line.png"),
  figure, device = "png", width = 16, height = 9, units = "cm", dpi = 200
)

########## Figura con Climatologia ############################################

A = matrix(
  c("left", "top", "center", "top", "right", "top", "left", "bot", "center", "bot", "right", "bot"),
  ncol = 2, nrow = 6, byrow = TRUE
)
  
  GG = vector("list", nq); names(GG) = quin
  
  for (q in 1:nq) {
    v = A[q,2]
    h = A[q,1]
    
    df = data.frame(
      lats = rep(rep(lats,nz), 4),
      plvl = rep(rep(plvl, each = ny), 4),
      uwind = c(
        array(DATA[,,q,c("2002", "2019")]), 
        array(apply(DATA[,,q,which(QBO[q,] == "QBO30-E")], c(1,2), mean)),
        array(apply(DATA[,,q,which(QBO[q,] == "QBO30-W")], c(1,2), mean))
      ),
      group = rep(c("2002","2019", "QBO30-E", "QBO30-W"), each = nz*ny)#,
      #colors = rep(QBO[q,], each = nz*ny)
    )
    
    GG[[q]] = ggplot() +
      geom_contour2(
        data = df, aes(x = lats, y = plvl, z = uwind, group = group, color = group),
        breaks = 0, size = 0.25
      ) +
      scale_color_manual(
        name = "Leyenda",
        values = c(
          "2019" = "blue", "2002"= "red",
          "QBO30-E" = "orange", "QBO30-W" = "purple"
        )
      ) +
      # scale_size_manual(
      #   values = c(
      #     "2019" = 0.2, "2002"= 0.2,
      #     "QBO30-E" = 0.1, "QBO30-W" = 0.1
      #   )
      # ) +
      # guides(
      #   size = "none"
      # ) +
      scale_y_level(
        name = NULL, breaks = c(100,70,50,30,20,10,7,5,3,2,1), sec.axis = dup_axis(name = NULL)
      ) +
      scale_x_latitude(
        breaks = seq(-60,0,20), limits = c(-70,10), sec.axis = dup_axis(name = NULL)
      ) +
      annotate("text", x = -60, y = 82, label = tit[q], color = "black", size = 2) +
      annotate("point", x = c(-60,0), y = c(10,30), size = 1) +
      theme_light() +
      theme(
        panel.ontop = TRUE, 
        panel.background = element_blank(), 
        panel.grid = element_line(size = 0.1, linetype = "dashed", color = "gray50"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        plot.title = element_blank(),
        legend.title = element_text(size = 5, face = "bold"),
        legend.text = element_text(size = 4),
        legend.position = "bottom",
        legend.margin = margin(c(0,0,0,0), unit = "cm"),
        legend.box.spacing = unit(0, units = "cm"),
        axis.title = element_blank(),
        # x-axis
        axis.text.x.bottom = tm.axis_text_x_bottom[[v]],
        axis.text.x.top = element_blank(),
        axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[v]],
        axis.ticks.x.top = tm.axis_ticks_x_top[[v]],
        # y-axis
        axis.text.y.left = tm.axis_text_y_left[[h]],
        axis.text.y.right = tm.axis_text_y_right[[h]],
        axis.ticks.y.left = tm.axis_ticks_y_left[[h]],
        axis.ticks.y.right = tm.axis_ticks_y_right[[h]],
        axis.text = element_text(size = 5)
      )
  }
  
  figure = ggarrange(plotlist = GG, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
  
  ggsave(
    here("outputs/ftropos/zero_wind_line_clim.png"),
    figure, device = "png", width = 16, height = 9, units = "cm", dpi = 200
  )









