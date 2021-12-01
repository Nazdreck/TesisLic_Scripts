########## Limpia el espacio de trabajo #######################################
rm(list=ls())
########## Paquetes y Funciones ###############################################
require(here)           # A simpler way to find your files
require(ncdf4)          # For manipulate .nc files
require(ggplot2)        # Data Visualisation
require(metR)           # Tools for Meteorological Analysis
require(colorspace)     # Paletas de Colores.
require(ggpubr)         # "ggplot2" Based Publication Ready Plots
require(gtools)

source(here("functions","QGEPflux.R"))   # Ask if year is leap

########## Predefiniendo variables ############################################
lons = seq(0, 357.5, 2.5); nx = length(lons)
lats = seq(0,-90,-2.5); ny = length(lats)
plvl = c(1000,700,500,300,200,100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
days = 1:30; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "1900-09-01")), 6,10)
years = c("2002", "2019"); na = length(years)
vars = c("t","u","v"); nv = length(vars)

DATA = array(
  double(nx*ny*nz*nt*na*nv), c(nx,ny,nz,nt,na,nv), 
  dimnames = list(lons = lons, lats = lats, plvl = plvl, dates = dates, years = years, vars = vars)
)

########## Extraccion de datos ################################################
for (v in 1:nv){
  for (a in years) {
    nc = nc_open(
      here("inputs/era5_plvls", paste0(vars[v], "_",a,".nc"))
    )
    DATA[,,,,a,v] = ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == lats[1]),
        1,
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "09-01")
      ), 
      count = c(-1,ny,-1,nt) 
    )
    nc_close(nc); rm(nc)
  }
  
  
}

# climatologia temperatura
nc = nc_open(
  here("inputs/era5_plvls/t_climatology.nc")
)
CLIM_t = ncvar_get(
  nc,
  start = c(
    1,
    which(ncvar_get(nc, "lats") == lats[1]),
    1,
    which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "08-01")
  ), 
  count = c(-1,ny,-1,nt) 
)
nc_close(nc); rm(nc)
dimnames(CLIM_t) = list(lons = lons, lats = lats, plvl = plvl, dates = dates)

print("Datos extraidos")

########## Flujo de Eliassen-Palm y Divergencia ###############################
FEP = array(
  double(ny*nz*nt*na*4), c(ny,nz,nt,na,4), 
  dimnames = list(lats = lats, plvl = plvl, dates = dates, years = years, vars = c("Fy", "Fp", "divF", "AcelTot"))
)

for (a in 1:na){ 
  for (t in 1:nt) {
    FEP[,,t,a,] = QGEPflux(dataset = DATA[,,,t,a,], longitudes = lons, latitudes = lats, pressure = plvl)
  }
}
FEP[c(1,ny),,,,] = NA

fdays = c("1-5","6-10","11-15","16-20", "21-25", "26-30"); nf = length(fdays)

F_proms = array(
  double(ny*nz*nf*na*3),c(ny,nz,nf,na,3),
  dimnames = list(lats = lats, plvl = plvl, dates = fdays, years = years, vars = c("Fy", "Fp", "AcelTot"))
)
ini = c(1,6,11,16,21,26)
for (i in 1:nf) {
  F_proms[,,i,,] = apply(FEP[,,ini[i]:(ini[i]+4),,c("Fy","Fp","AcelTot")], c(1,2,4,5), mean)
}

##### ESCALADO
### Taguchi y Hartman Scaling
S2 = array(rep((plvl/1000)^(1/2), each = ny), c(ny,nz))

sF = F_proms[,,,,c("Fy","Fp")]*0
for (a in 1:na) {
  for (t in 1:nf){
    sF[,,t,a,"Fy"] = F_proms[,,t,a,"Fy"]
    sF[,,t,a,"Fp"] = F_proms[,,t,a,"Fp"]/S2
  }
}


########## Anomalia de Temperatura y viento zonal medio #######################
Anom_T = DATA[,,,,,"t"]*0

for (a in 1:na) {
  Anom_T[,,,,a] = DATA[,,,,a,"t"] - CLIM_t
}
Anom_T = apply(Anom_T, 2:5, mean)

U = apply(DATA[,,,,,"u"], 2:5, mean)

T_proms = array(
  double(ny*nz*nf*na),c(ny,nz,nf,na),
  dimnames = list(lats = lats, plvl = plvl, dates = fdays, years = years)
)
U_proms = array(
  double(ny*nz*nf*na),c(ny,nz,nf,na),
  dimnames = list(lats = lats, plvl = plvl, dates = fdays, years = years)
)
for (i in 1:nf) {
  T_proms[,,i,] = apply(Anom_T[,,ini[i]:(ini[i]+4),], c(1,2,4), mean)
  U_proms[,,i,] = apply(U[,,ini[i]:(ini[i]+4),], c(1,2,4), mean)
}

########## Elementos para graficar ############################################
##### QUIEBRES y LABELS
brks = seq(-120,-15, 15); nb = length(brks)
labs = as.character(brks)
labs[which(even(1:nb))] = " "

##### COLORES
cols = sequential_hcl(nb-1, palette = "PuRd", rev = FALSE)

##### VECTORES MAXIMO Y MINIMO
# max_mag = max(Mag(FEP[,,,,"Fy"],FEP[,,,,"Fp"]), na.rm = TRUE)
# min_mag = 0.01*max_mag
max_mag = max(Mag(sF[,,,,"Fy"],sF[,,,,"Fp"]), na.rm = TRUE)
min_mag = 0.01*max_mag


##### Matriz iteradora
A = matrix(
  c(
    "1-5","6-10","11-15","16-20","21-25","26-30",
    "top","top","center","center","bot","bot",
    "left","right","left","right","left","right"
  ), 
  ncol = 3, nrow = 6
)

##### theming
# x axis
tm.axis_text_x_bottom = list("top" = element_text(), "center" = element_blank(), "bot" = element_blank())
tm.axis_text_x_top = list("top" = element_blank(), "center" = element_blank(), "bot" = element_text())
tm.axis_ticks_x_bottom = list("top" = element_line(), "center" = element_line(), "bot" = element_blank())
tm.axis_ticks_x_top = list("top" = element_blank(), "center" = element_line(), "bot" = element_line())
# y axis
tm.axis_text_y_left = list("left"  = element_text(), "right" = element_blank())
tm.axis_text_y_right = list("left"  = element_blank(), "right" = element_text())
tm.axis_ticks_y_left = list("left"  = element_line(), "right" = element_blank())
tm.axis_ticks_y_right = list("left"  = element_blank(), "right" = element_line())

########## Figuras ############################################################
for (a in 1:na) {
  GG = vector("list",4)
  for (t in 1:6) {
    tt = A[t,1]
    xx = A[t,2]
    yy = A[t,3]
    
    df = data.frame(
      lats = rep(lats, nz),
      plvl = rep(plvl, each = ny),
      Fy = array(sF[,,tt,a,"Fy"]),
      Fp = array(sF[,,tt,a,"Fp"]),
      divF = array(F_proms[,,tt,a,"AcelTot"]*86400),
      uwind = array(U_proms[,,tt,a]),
      temp = array(T_proms[,,tt,a])
    )
    
    GG[[t]] = ggplot(df) +
      aes(x = lats, y = plvl) +
      geom_contour_fill(aes(z = divF), breaks = brks) +
      scale_fill_gradientn(
        colors = cols, na.value = "white",
        breaks = brks,
        labels = labs,
        limits = range(brks)
      ) +
      guides(
        fill = guide_colorstrip(
          title = expression("d"~bar(u)~"/dt [m.s"^{-1}~".día"^{-1}~"]    "), barwidth = 7, barheight = 0.3,
          ticks = TRUE, ticks.color = "gray30", label.vjust = 2
        )
      ) +
      geom_contour2(aes(z = uwind), breaks = 0, linetype = "solid", size = .3, alpha = 0.6, color = "black") +
      geom_contour2(aes(z = uwind), breaks = seq(10,60,10), linetype = "solid", size = .15, alpha = 0.6, color = "gray20") +
      geom_contour2(aes(z = uwind), breaks = seq(-40,-10,10), linetype = "dashed", size = .15, alpha = 0.6, color = "gray20") +
      geom_text_contour(aes(z = uwind), breaks = seq(-40,60,10), size = 1.75, stroke = 0.05) +
      geom_contour2(aes(z = temp), breaks = seq(20,60,20), linetype = "solid", size = 0.2, color = "mediumblue") +
      #geom_contour2(aes(z = temp), breaks = -20, linetype = "dotdash", size = 0.2, color = "mediumblue") +
      geom_vector(
        aes(angle = atan2(Fp,Fy)*180/pi, mag = Mag(Fp,Fy)), color = "darkgreen", alpha = 0.9, pivot = 0, show.legend = FALSE,
        size = 0.4, arrow.length = unit(0.2, "cm"), min.mag = min_mag, direction = "cw", skip.x = 2, skip.y = 1
      ) +
      scale_mag(name = "F", max_size = 2.5, max = max_mag) +
      scale_x_latitude(name = NULL, breaks = seq(-15,-75,-15), sec.axis = dup_axis(name = NULL, breaks = seq(-15,-75,-15))) +
      scale_y_level(
        name = NULL, breaks = c(1000,700,500,300,200,100,70,50,30,20,10,7,5,3,2,1), sec.axis = dup_axis(name = NULL)
      ) +
      geom_vline(xintercept = c(-60,-30), size = 0.1, color = "gray30", linetype = "dashed") +
      geom_hline(yintercept = c(100,10), size = 0.1, color = "gray30", linetype = "dashed") +
      annotate(geom = "text", x = -83, y = 1.4, label = fdays[t], size = 2.5) +
      theme_light() +
      theme(
        # panel.ontop = TRUE, 
        # panel.background = element_blank(), 
        panel.grid = element_blank(), 
        # panel.grid.minor = element_blank(),
        plot.margin = margin(c(0,0,0,0), unit = "cm"),
        # plot.title = element_text(size = 6, face = "bold", vjust = -8, hjust = 0.03, margin = margin(c(-0.19,0,0,0), unit = "cm")),
        plot.title = element_blank(),
        legend.title = element_text(size = 5, face = "bold"), 
        legend.text = element_text(size = 4), 
        legend.position = "bottom",
        legend.margin = margin(c(0,0,0,0), unit = "cm"),
        legend.box.spacing = unit(0, units = "cm"),
        axis.title = element_blank(),
        # x-axis
        axis.text.x.bottom = tm.axis_text_x_bottom[[xx]],
        axis.text.x.top = tm.axis_text_x_top[[xx]],
        axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[xx]],
        axis.ticks.x.top = tm.axis_ticks_x_top[[xx]],
        # y-axis
        axis.text.y.left = tm.axis_text_y_left[[yy]],
        axis.text.y.right = tm.axis_text_y_right[[yy]],
        axis.ticks.y.left = tm.axis_ticks_y_left[[yy]],
        axis.ticks.y.right = tm.axis_ticks_y_right[[yy]],
        axis.text = element_text(size = 5, color = "gray50")
      )
    

  }
  
  figure = ggarrange(plotlist = GG, ncol = 2, nrow = 3, common.legend = TRUE, legend = "bottom")
  
  
  ggsave(
    here("outputs/ftropos", paste0("FEP_plvl-lats_", years[a],".png")),
    figure, device = "png", width = 12, height = 18, units = "cm", dpi = 200
  )
}

