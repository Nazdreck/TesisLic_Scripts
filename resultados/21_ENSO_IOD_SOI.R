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

########## Extrayendo datos ###################################################
NINO34 = read.table(
  here("inputs/ENSO/NINO34"), sep = "", dec = ".", row.names = 1
)

SOI = read.table(
  here("inputs/ENSO/SOI"), sep = "", dec = ".", row.names = 1
)

DOI =  read.table(
  here("inputs/ENSO/DOI"), sep = "", dec = ".", row.names = 1
)

DOI = array(t(as.matrix(DOI)))
NINO34 = NINO34[[1]]
SOI = SOI[[1]]
########## Graficos ###########################################################
# NINO34

df1 = data.frame(
  time = 1:372,
  var = NINO34
)

ypos = NINO34
ypos[which(ypos <= 0)] = 0
yneg = NINO34
yneg[which(yneg >= 0)] = 0
rib1 = list(ypos,yneg)

GG1 = ggplot(df1) +
  aes(x = time, y = var) +
  geom_ribbon(aes(ymin = 0, ymax = rib1[[1]]), fill = "red") +
  geom_ribbon(aes(ymin = rib1[[2]], ymax = 0), fill = "blue") +
  geom_line(size = 0.35, alpha = 0.7) + 
  annotate("rect", xmin = 354, xmax = 357, ymin = -3, ymax = 3, alpha = 0.4, fill = "green") +
  annotate("rect", xmin = 151, xmax = 154, ymin = -3, ymax = 3, alpha = 0.4, fill = "green") +
  scale_y_continuous(name = "Anomalía TSM [°C]", breaks = -3:3, limits = c(-3,3), expand = c(0,0)) +
  scale_x_continuous(name = NULL, breaks = c(seq(1,360,12*5),360), labels = seq(1990,2020,5), limits = range(109:372), expand = c(0,0)) +
  labs(title = "El Niño 3.4") +
  theme_light() +
  theme(
    plot.margin = margin(t = 0, r = 0.3, b = 0, l = 0.1, unit = "cm"),
    plot.title = element_text(face = "bold", size = 6, vjust = -3),
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 6)
  )

# SOI

df2 = data.frame(
  time = 1:372,
  var = SOI
)

ypos = SOI
ypos[which(ypos <= 0)] = 0
yneg = SOI
yneg[which(yneg >= 0)] = 0
rib2 = list(ypos,yneg)

GG2 = ggplot(df2) +
  aes(x = time, y = var) +
  geom_ribbon(aes(ymin = 0, ymax = rib2[[1]]), fill = "blue") +
  geom_ribbon(aes(ymin = rib2[[2]], ymax = 0), fill = "red") +
  geom_line(size = 0.35, alpha = 0.7) + 
  annotate("rect", xmin = 354, xmax = 357, ymin = -6, ymax = 6, alpha = 0.4, fill = "green") +
  annotate("rect", xmin = 151, xmax = 154, ymin = -6, ymax = 6, alpha = 0.4, fill = "green") +
  scale_y_continuous(name = NULL, breaks = seq(-6,6,2), limits = c(-6,6), expand = c(0,0)) +
  scale_x_continuous(name = NULL, breaks = c(seq(1,360,12*5),360), labels = seq(1990,2020,5), limits = range(109:372), expand = c(0,0)) +
  labs(title = "Índice Oscilacion del Sur (SOI)") +
  theme_light() +
  theme(
    plot.margin = margin(t = 0, r = 0.3, b = 0, l = 0.1, unit = "cm"),
    plot.title = element_text(face = "bold", size = 6, vjust = -3),
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 6)
  )
  
# DOI

df3 = data.frame(
  time = 1:360,
  var = DOI
)

ypos = DOI
ypos[which(ypos <= 0)] = 0
yneg = DOI
yneg[which(yneg >= 0)] = 0
rib3 = list(ypos,yneg)

GG3 = ggplot(df3) +
  aes(x = time, y = var) +
  geom_ribbon(aes(ymin = 0, ymax = rib3[[1]]), fill = "red") +
  geom_ribbon(aes(ymin = rib3[[2]], ymax = 0), fill = "blue") +
  geom_line(size = 0.35, alpha = 0.7) + 
  annotate("rect", xmin = 354, xmax = 357, ymin = -1.5, ymax = 1.5, alpha = 0.4, fill = "green") +
  annotate("rect", xmin = 151, xmax = 154, ymin = -1.5, ymax = 1.5, alpha = 0.4, fill = "green") +
  scale_y_continuous(name = "Anomalía TSM [°C]", breaks = seq(-1.5,1.5,0.5), limits = c(-1.5,1.5), expand = c(0,0)) +
  scale_x_continuous(name = NULL, breaks = c(seq(1,360,12*5),360), labels = seq(1990,2020,5), limits = range(121:360), expand = c(0,0)) +
  labs(title = "Dipolo del Oceano Índico (IOD)") +
  theme_light() +
  theme(
    plot.margin = margin(t = 0, r = 0.3, b = 0, l = 0.1, unit = "cm"),
    plot.title = element_text(face = "bold", size = 6, vjust = -3),
    axis.text = element_text(size = 5),
    axis.title.x = element_blank(),
    axis.title = element_text(size = 6)
  )  
  
figure = ggarrange(
  plotlist = list(GG1,GG2,GG3), nrow = 3, ncol = 1, align = "v"
)

ggsave(
  here("outputs/ftropos/TSM_indexes.png"),
  figure, device = "png", width = 12, height = 9, units = "cm", dpi = 200
)

  
  
  