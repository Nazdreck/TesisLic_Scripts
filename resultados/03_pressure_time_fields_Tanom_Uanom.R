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
vars = c("t","u"); nv = length(vars)
days = 1:365; nt = length(days)
dates = substr(as.character(as.Date(days-1, origin = "2000-01-01")), 6, 10)
plvl = c(1000,700,500,300,200,100,70,50,30,20,10,7,5,3,2,1); nz = length(plvl)
stats = c("2002", "2019", "climatology"); ns = length(stats)

months = c("01","02","03","04","05","06","07","08","09","10","11","12"); nm = length(months)
years = c(2002,2019); na = length(years)

########## Extrayendo datos ###################################################
DATA = array(double(nt*nz), c(nt,nz,ns,nv), dimnames = list(dates = dates, plvl = plvl, stats = stats, vars = vars))

for (v in vars) {
  if (v == "t") lt = -1 else lt = 1
  for (s in stats) {
    nc = nc_open(
      here("inputs/era5_plvls", paste0(v, "_", s, ".nc"))
    )
    
    data_extracted = ncvar_get(
      nc,
      start = c(
        1,
        which(ncvar_get(nc, "lats") == -60),
        1, 1
      ),
      count = c(
        c(-1,lt,-1,-1)
      )
    )
    nc_close(nc); rm(nc)
    
    if (v == "t") data_extracted = t(apply(data_extracted, c(3,4), mean))
    if (v == "u") data_extracted = t(apply(data_extracted, c(2,3), mean))
    
    DATA[,,s,v] = data_extracted
  }
}
rm(data_extracted)

# Carga los dias con minimos de flujo de calor
load("inputs/heatflux_mindays.RData")

########## Anomalias ##########################################################
ANOM = array(double(nt*nz*na*nv), c(nt,nz,na,nv), list(dates = dates, plvl = plvl, years = years, vars = vars))

for (v in vars) {
  for (a in as.character(years)) {
    ANOM[,,a,v] = DATA[,,a,v] - DATA[,,"climatology",v]
  }
}

########## Elementos para graficar ############################################
mini = c(1,32,60,91,121,152,182,213,244,274,305,335,365)
mname = c("Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic", "Ene")

brks = list(t = seq(-40,40,5), u = seq(-100,100,10)); nb = c(t = length(brks$t), u = length(brks$u))
labs = list(t = as.character(brks$t), u = as.character(brks$u))
labs$t[which(even(1:nb["t"]))] = " "
labs$u[which(even(1:nb["u"]))] = " "

colors = list(
  t = c(
    sequential_hcl((nb["t"]-1)/2, palette = "Blues", rev = FALSE), # YlGnBu / Blues /
    sequential_hcl((nb["t"]-1)/2, palette = "Reds", rev = TRUE) # YlOrRd / Reds / OrRd
  ),
  u = c(
    sequential_hcl((nb["u"]-1)/2, palette = "PuRd", rev = FALSE), # PuRd / Purples / RdPu / BuPu
    sequential_hcl((nb["u"]-1)/2, palette = "Greens", rev = TRUE) # BuGn / Greens /
  )
)

legend_title = c(t = "K   ", u = "m/s   ")
plot_title = list(t = expression(bar(T)~ "(60-90°S)"), u = expression(bar(u)~ "(60°S)"))
hfmin = list("2002" = mins2002, "2019" = mins2019)

# Theming
tm.plot_margin = list(
  "t" = list(
    "2002" = margin(t = 0, r = 0.075, b = -0.01, l = 0.075, unit = "cm"), 
    "2019" = margin(t = 0.07, r = 0.075, b = 0.12, l = 0.075, unit = "cm")
  ),
  "u" = list(
    "2002" = margin(t = 0, r = 0.075, b = -0.01, l = 0.075, unit = "cm"), 
    "2019" = margin(t = 0.07, r = 0.075, b = 0.12, l = 0.075, unit = "cm")
  )
)
tm.plot_title = list("2002" = element_text(size = 7, face = "bold", vjust = -1.5, margin = margin(c(-0.1,0,0,0), unit = "cm")), "2019" = element_blank())
tm.axis_text_x_bottom = list("2002" = element_text(), "2019" = element_blank())
tm.axis_ticks_x_bottom = list("2002" = element_line(), "2019" = element_blank())
tm.axis_ticks_x_top = list("2002" = element_blank(), "2019" = element_line())
tm.axis_text_y_left = list("t" = element_text(), "u" = element_blank())
tm.axis_text_y_right = list("t" = element_blank(), "u" = element_text())
tm.axis_title_y_left = list("t" = element_text(margin = margin(r = -0.1, unit = "cm")), "u" = element_blank())
tm.axis_title_y_right = list("t" = element_blank(), "u" = element_text(margin = margin(l = -0.1, unit = "cm")))
tm.axis_ticks_y_left = list("t" = element_line(), "u" = element_blank())
tm.axis_ticks_y_right = list("t" = element_blank(), "u" = element_line())
########## Graficos ###########################################################
GG = vector("list", nv); names(GG) = vars

for (v in vars) {
  gg_a = vector("list",na); names(gg_a) = years
  for (a in as.character(years)) {
    
    df = data.frame(
      days = rep(days, nz),
      plvl = rep(plvl, each = nt),
      var = array(ANOM[,,a,v])
    )
    
    gg_a[[a]] = ggplot(df) +
      aes(x = days, y = plvl, z = var) +
      geom_contour_fill(breaks = brks[[v]]) +
      scale_fill_gradientn(
        colors = colors[[v]],
        breaks = brks[[v]],
        labels = labs[[v]],
        limits = range(brks[[v]])
      ) +
      guides(
        fill = guide_colorstrip(
          title = legend_title[v], ticks = TRUE, ticks.color = "gray30", barwidth = 10, barheight = 0.3,
          label.vjust = 2
        )
      ) +
      geom_contour2(breaks = 0, size = 0.1, color = "gray40", alpha = 0.6) +
      scale_x_continuous(
        name = NULL, limits = range(121:335), breaks = mini, labels = mname, expand = c(0,0), sec.axis = dup_axis(name = NULL)
      ) +
      scale_y_level(
        name = "hPa", breaks = c(1000,500,200,100,50,20,10,5,2,1), sec.axis = dup_axis(name = "hPa")
      ) +
      geom_vline(xintercept = hfmin[[a]], color = "blue", linetype = "dashed", size = 0.1) +
      geom_hline(yintercept = c(100,10), color = "gray50", linetype = "dashed", size = 0.05) +
      annotate("text", x = 130, y = 700, label = a, color = "black", size = 1.75) +
      labs(title = plot_title[[v]]) +
      theme_light() +
      theme(
        plot.title = tm.plot_title[[a]],
        plot.margin = tm.plot_margin[[v]][[a]],
        legend.title = element_text(size = 5, face = "bold"),
        legend.text = element_text(size = 4),
        legend.position = "bottom",
        legend.margin = margin(c(0,0,0,0), unit = "cm"),
        legend.box.spacing = unit(0, units = "cm"),
        axis.title = element_text(size = 5, face = "bold", color = "gray50"),
        axis.text.x.bottom = tm.axis_text_x_bottom[[a]],
        axis.text.x.top = element_blank(),
        axis.ticks.x.bottom = tm.axis_ticks_x_bottom[[a]],
        axis.ticks.x.top = tm.axis_ticks_x_top[[a]],
        axis.title.y.left = tm.axis_title_y_left[[v]],
        axis.title.y.right = tm.axis_title_y_right[[v]],
        axis.text.y.left = tm.axis_text_y_left[[v]],
        axis.text.y.right = tm.axis_text_y_right[[v]],
        axis.ticks.y.left = tm.axis_ticks_y_left[[v]],
        axis.ticks.y.right = tm.axis_ticks_y_right[[v]],
        axis.text = element_text(size = 5)

      )
  }
  
  GG[[v]] = ggarrange(plotlist = gg_a, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
}

figure = ggarrange(plotlist = GG, ncol = 2, nrow = 1)

ggsave(
  here("outputs/ftropos/T-U_vertical_anom.png"),
  figure, device = "png", width = 16, height = 9, units = "cm", dpi = 200
)
