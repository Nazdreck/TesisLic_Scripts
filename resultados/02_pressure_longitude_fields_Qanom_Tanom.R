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
lons = seq(0,357.5,2.5); nx = length(lons)
lats = seq(-60,-90,-2.5); ny = length(lats)
plvl = c(100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
days = 1:61; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "2000-09-01")), 6, 10)
months = c("09","10"); nm = length(months)
stats = c("2002", "2019", "climatology"); ns = length(stats)
years = c(2002,2019); na = length(years)
vars = c("t","pv"); nv = length(vars)

########## Extrayendo datos ###################################################
DATA = array(double(nx*ny*nz*nt*ns*nv), c(nx,ny,nz,nt,ns,nv), dimnames = list(
  lons = lons, lats = lats, plvl = plvl, dates = dates, stats = stats, vars = vars)
)

for (v in vars) {
  for (s in stats) {
    nc = nc_open(
      here("inputs/era5_plvls", paste0(v, "_", s, ".nc"))
    )
    
    data_extracted = ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == lats[1]),
        which(ncvar_get(nc, "plvl") == plvl[1]),
        which(substr(as.character(as.Date(ncvar_get(nc, "time")-1, origin = "1900-01-01")),6,10) == "09-01")
      ),
      count = c(
        c(-1,-1,-1,nt)
      )
    )
    nc_close(nc); rm(nc)
    
    DATA[,,,,s,v] = data_extracted
    
  }
}
rm(data_extracted)
########## Preparando datos ###################################################
TEMP = DATA[,,,,c("2002","2019"),"t"]
anom_PTVT = DATA[,,,,c("2002","2019"),"pv"]
norm_PTVT = DATA[,,,,c("2002","2019"),"pv"]

for (a in as.character(years)) {
  TEMP[,,,,a] = TEMP[,,,,a] - DATA[,,,,"climatology","t"]
  anom_PTVT[,,,,a] = anom_PTVT[,,,,a] - DATA[,,,,"climatology","pv"]
  norm_PTVT[,,,,a] = (norm_PTVT[,,,,a] - DATA[,,,,"climatology","pv"])/DATA[,,,,"climatology","pv"]
}

TEMP = apply(TEMP, c(1,3,4,5), mean)
anom_PTVT = apply(anom_PTVT, c(1,3,4,5), mean)
norm_PTVT = apply(norm_PTVT, c(1,3,4,5), mean)
PTVT = apply(DATA[,,,,c("2002","2019"),"pv"], c(1,3,4,5), mean)


dts = c("-14d","0d","+14d"); nd = length(dts)
mm = list("2002" = c(-14,0,14) + 26, "2019" = c(-14,0,14) + 20)
dd = list(
  "2002" = matrix(vector("character",15), nrow = 5, ncol = 3),
  "2019" = matrix(vector("character",15), nrow = 5, ncol = 3)
)
for (a in as.character(years)) {
  for (i in 1:3) {
    dd[[a]][,i] = dates[(mm[[a]][i]-2):(mm[[a]][i]+2)]
  }
  
}

dTEMP = array(double(nx*nz*nd*na), c(nx,nz,nd,na), dimnames = list(lons = lons, plvl = plvl, days = dts, years = years))
dPTVT = array(double(nx*nz*nd*na), c(nx,nz,nd,na), dimnames = list(lons = lons, plvl = plvl, days = dts, years = years))
danom_PTVT = array(double(nx*nz*nd*na), c(nx,nz,nd,na), dimnames = list(lons = lons, plvl = plvl, days = dts, years = years))
dnorm_PTVT = array(double(nx*nz*nd*na), c(nx,nz,nd,na), dimnames = list(lons = lons, plvl = plvl, days = dts, years = years))

for (a in as.character(years)) {
  for (d in 1:nd) {
    dTEMP[,,d,a] = apply(TEMP[,,dd[[a]][,d],a], c(1,2), mean)
    dPTVT[,,d,a] = apply(PTVT[,,dd[[a]][,d],a], c(1,2), mean)
    danom_PTVT[,,d,a] = apply(anom_PTVT[,,dd[[a]][,d],a], c(1,2), mean)
    dnorm_PTVT[,,d,a] = apply(norm_PTVT[,,dd[[a]][,d],a], c(1,2), mean)
  }
}

dPTVT = dPTVT/10^-6
danom_PTVT = danom_PTVT/10^-6


# m02 = which(dates == "09-26"); m19 = which(dates == "09-20")
# dd = list("2002" = dates[(m02-2):(m02+2)], "2019" = dates[(m19-2):(m19+2)])
# 
# dTEMP = array(double(nx*nz*na), c(nx,nz,na), dimnames = list(lons = lons, plvl = plvl, years = years))
# dPTVT = array(double(nx*nz*na), c(nx,nz,na), dimnames = list(lons = lons, plvl = plvl, years = years))
# danom_PTVT = array(double(nx*nz*na), c(nx,nz,na), dimnames = list(lons = lons, plvl = plvl, years = years))
# 
# for (a in as.character(years)) {
#   dTEMP[,,a] = apply(TEMP[,,dd[[a]],a], c(1,2), mean)
#   dPTVT[,,a] = apply(PTVT[,,dd[[a]],a], c(1,2), mean)
#   danom_PTVT[,,a] = apply(anom_PTVT[,,dd[[a]],a], c(1,2), mean)
# }
# dPTVT = dPTVT/10^-6
# danom_PTVT = danom_PTVT/10^-6
########## Elementos para graficar ############################################
brks_t = seq(-50,50,5); nb = length(brks_t)
brks_pv = c(-10000,-7000,-5000,-3000,-2000,-1000,-700,-500,-300,-200,-100,-70,-50,-30,-20,-10)
brks_pv_p = c(10,20,30,50,70,100,200,300,500,700,1000,2000,3000,5000)
brks_pv_n = c(-1000,-700,-500,-300,-200,-100,-70,-50,-30,-20,-10)

labs = as.character(brks_t)
labs[which(even(1:nb))] = " "

cols = c(
  sequential_hcl((nb-1)/2, palette = "Blues", rev = FALSE), # YlGnBu / Blues /
  sequential_hcl((nb-1)/2, palette = "Reds", rev = TRUE) # YlOrRd / Reds / OrRd
)

A = matrix(
  c(rep(2002,3), rep(2019,3), rep(1:3, 2)), ncol = 2, nrow = 6
)

##### theming
# x axis
tm.axis_text_x_bottom = list("2002" = element_text(), "2019" = element_blank())
tm.axis_ticks_x_bottom = list("2002" = element_line(), "2019" = element_blank())
tm.axis_ticks_x_top = list("2002" = element_blank(), "2019" = element_line())
# y axis
#tm.axis_title_y_left = list("t" = element_text(margin = margin(r = -0.1, unit = "cm")), "u" = element_blank())
#tm.axis_title_y_right = list("t" = element_blank(), "u" = element_text(margin = margin(l = -0.1, unit = "cm")))
tm.axis_text_y_left = list("1" = element_text(), "2" = element_blank(), "3" = element_blank())
tm.axis_text_y_right = list("1" = element_blank(), "2" = element_blank(), "3" = element_text())
tm.axis_ticks_y_left = list("2" = element_line(), "2" = element_blank(), "3" = element_blank())
tm.axis_ticks_y_right = list("2" = element_blank(), "2" = element_blank(), "3" = element_line())
########## Graficos ###########################################################
GG = vector("list", na*nd)
for (i in 1:(na*nd)) {
  a = as.character(A[i,1])
  d = A[i,2]
  
  df = data.frame(
    lons = rep(ConvertLongitude(c(lons[73:144],lons[1:72])), nz),
    plvl = rep(plvl, each = nx),
    temp = array(abind(dTEMP[73:144,,d,a], dTEMP[1:72,,d,a], along = 1)),
    #pv = array(abind(dPTVT[73:144,,d,a], dPTVT[1:72,,d,a], along = 1 ))
    #pv = array(abind(danom_PTVT[73:144,,d,a], danom_PTVT[1:72,,d,a], along = 1 ))
    pv = array(abind(dnorm_PTVT[73:144,,d,a], dnorm_PTVT[1:72,,d,a], along = 1 ))
    
  )
  
  GG[[i]] = ggplot(df) +
    aes(x = lons, y = plvl) +
    geom_contour_fill(aes(z = temp), breaks = brks_t) +
    scale_fill_gradientn(
      colors = cols,
      breaks = brks_t,
      labels = labs,
      limits = range(brks_t)
    ) +
    guides(
      fill = guide_colorstrip(
        title = "T [K]   ", ticks = TRUE, ticks.color = "gray30", barwidth = 10, barheight = 0.3, label.vjust = 2
      )
    ) +
    # geom_contour2(aes(z = pv), breaks = brks_pv_p, linetype = "solid", size = 0.1) +
    # geom_contour2(aes(z = pv), breaks = 0, linetype = "solid", color = "gray80", size = 0.2) +
    # geom_contour2(aes(z = pv), breaks = brks_pv_n, linetype = "dashed", size = 0.1) +
    # geom_text_contour(aes(z = pv), breaks = c(brks_pv_n,0,brks_pv_p), stroke = 0.1, size = 2) +
    ###
    # geom_contour2(aes(z = pv), breaks = brks_pv, size = 0.1) +
    # geom_text_contour(aes(z = pv), breaks = brks_pv, stroke = 0.1, size = 1.75) +
    ###
    geom_contour2(aes(z = pv), breaks = seq(0.2,1,0.2), linetype = "solid", size = 0.1) +
    geom_contour2(aes(z = pv), breaks = 0, linetype = "solid", color = "gray80", size = 0.2) +
    geom_contour2(aes(z = pv), breaks = seq(-1,-0.2,0.2), linetype = "dashed", size = 0.1) +
    geom_text_contour(aes(z = pv), breaks = seq(-1,1,0.2), stroke = 0.1, size = 2.2) +
    scale_y_level(
      name = NULL, breaks = c(100,70,50,30,20,10,7,5,3,2,1), sec.axis = dup_axis(name = NULL)
    ) +
    scale_x_longitude(
      breaks = c(-120,-60,0,60,120), sec.axis = dup_axis(name = NULL)
    ) +
    annotate("text", x = -130, y = 80, label = paste0("(",a,") ", "[",dts[d],"]" ), color = "black", size = 1.75) +
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
      #axis.title = element_text(size = 5, face = "bold", color = "gray50"),
      axis.text.x.bottom = tm.axis_text_x_bottom[[a]],
      axis.text.x.top = element_blank(),
      axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[a]],
      axis.ticks.x.top = tm.axis_ticks_x_top[[a]],
      # axis.title.y.left = tm.axis_title_y_left[[v]],
      # axis.title.y.right = tm.axis_title_y_right[[v]],
      axis.text.y.left = tm.axis_text_y_left[[d]],
      axis.text.y.right = tm.axis_text_y_right[[d]],
      axis.ticks.y.left = tm.axis_ticks_y_left[[d]],
      axis.ticks.y.right = tm.axis_ticks_y_right[[d]],
      axis.text = element_text(size = 5),
    )
} 

figure = ggarrange(plotlist = GG, ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

ggsave(
  here("outputs/ftropos/T-PV_zx_norm.png"),
  figure, device = "png", width = 16, height = 9, units = "cm", dpi = 200
)














