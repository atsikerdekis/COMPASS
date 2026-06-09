#################
### LIBRARIES ###
#################
library(ggplot2)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)

#############
### INPUT ###
#############
expname1 <- "it41"
exptype1 <- "AER"
expname2 <- "b2th"
exptype2 <- "HAM"
expnames <- c(expname1, expname2)
sDate   <- as.Date("2018-12-15")
eDate   <- as.Date("2018-12-31")
seqDate <- seq.Date(sDate, eDate, "day")

#################
### READ DATA ###
#################
mydata <- NULL
for (e in 1:length(expnames)) {
  for (d in 1:length(seqDate)) { 
    temp <- read.table(paste0("cams2_35_bis/data/massdia/",expnames[e],"/massdia_chem__",expnames[e],"_",gsub("-","",seqDate[d]),"00.txt"), head=T)
    temp$date <- seqDate[d]
    temp$datetime <- as.POSIXct(paste0(temp$date," ",temp$SIM_HOUR,":00:00"), tz="GMT")
    temp$SIM_HOUR <- NULL
    if (!is.null(temp$SEDM_FLX)) { temp$DDEP_FLX <- temp$DDEP_FLX + temp$SEDM_FLX }
    temp$SEDM_FLX <- NULL
    temp$expname  <- expnames[e]
    mydata <- rbind(mydata,temp)
  }
}

####################
### POST-PROCESS ###
####################
### Select only the 24th timestep, which includes the tendency for all previous timesteps in that day!!!
### See discussion in slack with Vincent... well it will probably disappear in 3 months from now 28/08/2025
### Units are in Tg/day in the SIM_HOUR==24
### Units are in Tg/6hours for the first 6 hours of the day in SIM_HOUR==6
mydata <- mydata[which(format(mydata$datetime,"%H")=="00"),]

# convert to data.table
mydata <- setDT(mydata)

#"Sea-salt_1","Sea-salt_2","Sea-salt_3",
#"Desert-dust_1","Desert-dust_2","Desert-dust_3",
#"Organic-matter_A","Organic-matter_B",
#"Black-carbon_A","Black-carbon_B",
#"Sulphate_SO4",
#"Nitrate_1","Nitrate_2",
#"Ammonium",
#"SecOrg_Bio","SecOrg_Anth"

#"POM_KI","BC_KI","SOA_KI",
#"DU_AI",
#"DU_CI",
#"SO4_NS","SOA_NS",
#"POM_KS","BC_KS","SO4_KS","SOA_KS",
#"DU_AS","SS_AS","POM_AS","BC_AS","SO4_AS","AM_AS","NI_AS","SOA_AS",
#"DU_CS","SS_CS","POM_CS","BC_CS","SO4_CS"

# define groups: NAME -> group
group_map <- list(
  DU_AER = c("Desert-dust_1","Desert-dust_2","Desert-dust_3"),
  SS_AER = c("Sea-salt_1","Sea-salt_2","Sea-salt_3"),
  OM_AER = c("Organic-matter_A","Organic-matter_B"),
  BC_AER = c("Black-carbon_A","Black-carbon_B"),
  SU_AER = c("Sulphate_SO4"),
  NI_AER = c("Nitrate_1", "Nitrate_2"),
  AM_AER = c("Ammonium"),
  SOA_AER = c("SecOrg_Bio","SecOrg_Anth"),
  DU_HAM = c("DU_AI","DU_CI","DU_AS","DU_CS"),
  SS_HAM = c("SS_AS","SS_CS"),
  OM_HAM = c("POM_KI","POM_KS","POM_AS","POM_CS"),
  BC_HAM = c("BC_KI","BC_KS","BC_AS","BC_CS"),
  SU_HAM = c("SO4_NS","SO4_KS","SO4_AS","SO4_CS"),
  NI_HAM = c("NI_AS"),
  AM_HAM = c("AM_AS"),
  SOA_HAM = c("SOA_KI","SOA_NS","SOA_KS","SOA_AS")
)

# function to compute grouped means
aggregate_groups <- function(dt, group_map) {
  out_list <- lapply(names(group_map), function(newname) {
    members <- group_map[[newname]]
    
    tmp <- dt[NAME %in% members, 
              lapply(.SD, mean, na.rm = TRUE), 
              by = .(expname, datetime),
              .SDcols = c("TOT_MASS","TRP_MASS","NEG_MASS","EMIS_FLX",
                          "DDEP_FLX","WDEP_FLX","CHEM_TND","CHEM_TRO",
                          "NEGA_FIX","FLUX_ERR","EMIS_ATM",
                          "MASS_CHG","RESIDUAL","MASS_CHG_r","RESIDUAL_r")]
    tmp[, NAME := newname]
    return(tmp)
  })
  rbindlist(out_list)
}

# run aggregation
df_grouped <- aggregate_groups(mydata, group_map)
df_grouped <- df_grouped[,c("NAME","datetime","CHEM_TND","DDEP_FLX","EMIS_FLX","WDEP_FLX","NEGA_FIX","TOT_MASS")]

#> mydata_ham
#        NAME       date   datetime      CHEM_TND      DDEP_FLX     EMIS_FLX      WDEP_FLX     NEGA_FIX
#541      SS1 2018-12-01 2018-12-02 -0.0004079739 -2.547367e-01 9.551235e-01 -1.240406e+00 1.513907e-03
#542      SS2 2018-12-01 2018-12-02  0.0000000000 -1.085939e+01 2.528663e+01 -1.531483e+01 2.118855e-02
#543      SS3 2018-12-01 2018-12-02  0.0000000000 -5.654060e+01 9.324447e+01 -1.677105e+01 2.888850e-02


df <- as.data.frame(df_grouped)

df_long <- df %>%
  pivot_longer(
    cols = where(is.numeric),   # all numeric columns
    names_to = "variable",
    values_to = "value"
  )

df_means <- df_long %>%
  group_by(NAME, variable) %>%
  summarise(period_mean = mean(value, na.rm = TRUE), .groups = "drop")

### UNITS from Tg/day to Tg/year
df_means$period_mean <- df_means$period_mean*365

# Ensure NAME has a defined order
df_means$NAME <- factor(df_means$NAME, levels = unique(df_means$NAME))



# --- prepare data ----------------------------------------------------------
lvls <- c("DU","SS","OM","BC","SU","NI","AM","SOA")  # check 'SU' vs your actual label
df_plot <- df_means %>%
  mutate(
    TYPE = ifelse(stringr::str_detect(NAME, "HAM"), "HAM", "AER"),
    SPEC = stringr::str_remove(NAME, "_.*$") |> factor(levels = lvls)
  )


# --- stacked (exclude TOT_MASS) and rename variables -----------------------
x_levels <- c("HAM", "AER", "SPACER", "HAM_TOT", "AER_TOT")
x_labels <- setNames(
  c("Processes (HAM)", "Processes (AER)", "", "Total Mass (HAM)", "Total Mass (AER)"),
  x_levels
)

df_stacked <- df_plot %>%
  filter(variable != "TOT_MASS") %>%
  mutate(variable = recode(variable,
                           "CHEM_TND" = "CHEM",
                           "DDEP_FLX"  = "DDEP",
                           "WDEP_FLX"  = "WDEP",
                           "EMIS_FLX"  = "EMIS",
                           "NEGA_FIX"  = "NEG_FIX"),
         x_slot = factor(TYPE, levels = x_levels)
  )

# --- sum of processes (one row per NAME/TYPE/SPEC) -------------------------
df_sums <- df_stacked %>%
  group_by(NAME, TYPE, SPEC, x_slot) %>%
  summarise(total_proc = sum(period_mean), .groups = "drop")

# --- TOT_MASS bars ---------------------------------------------------------
df_tot <- df_plot %>%
  filter(variable == "TOT_MASS") %>%
  mutate(x_slot = factor(paste0(TYPE, "_TOT"), levels = x_levels))

# --- palette ----------------------------------------------------------------
my_colors <- c(
  "WDEP"    = "#4A90E2",
  "DDEP"    = "#593516",
  "EMIS"    = "#8E44AD",
  "CHEM"    = "#E74C3C",
  "NEG_FIX" = "#109c0b"
)

# --- plot -------------------------------------------------------------------
# --- calculate symmetric ranges per SPEC ---
df_range <- df_plot %>%
  group_by(SPEC) %>%
  summarise(max_val = max(abs(period_mean), na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(max_val))  # drop specs with only NAs

# --- make sure x_slot factors are aligned everywhere ---
df_stacked$x_slot <- factor(df_stacked$x_slot, levels = x_levels)
df_tot$x_slot     <- factor(df_tot$x_slot,     levels = x_levels)
df_sums$x_slot    <- factor(df_sums$x_slot,    levels = x_levels)

# --- plot ---
p <- ggplot() +
  geom_bar(
    data = df_stacked,
    aes(x = x_slot, y = period_mean, fill = variable),
    stat = "identity", position = "stack", width = 0.9
  ) +
  geom_crossbar(
    data = df_sums,
    aes(x = x_slot, y = total_proc, ymin = total_proc, ymax = total_proc),
    width = 0.6, color = "black", size = 0.7
  ) +
  geom_text(
    data = df_sums,
    aes(x = x_slot, y = total_proc, label = round(total_proc, 1)),
    vjust = -0.7, size = 3.2, fontface = "bold"
  ) +
  geom_bar(
    data = df_tot,
    aes(x = x_slot, y = period_mean),
    stat = "identity", fill = "grey70", width = 0.9
  ) +
  geom_text(
    data = df_tot,
    aes(x = x_slot, y = period_mean, label = round(period_mean, 1)),
    vjust = -0.7, size = 3.2, fontface = "bold"
  ) +
  facet_wrap(~SPEC, ncol = 4, scales = "free_y") +
  scale_x_discrete(
    limits = x_levels,
    labels = c("HAM", "AER", "", "HAM", "AER"),  # blank label = gap, no tick
    drop = FALSE,
    expand = expansion(mult = c(0, 0.05))       # creates a small spacing
  ) +
  scale_fill_manual(values = my_colors, name = "Gain/Loss Process") +
  labs(x = NULL, 
       y = "Fluxes (Tg/year)",
       title = paste0("Global diagnostics of ",expname1," (",exptype1,") and ",expname2," (",exptype2,") for ",sDate," to ",eDate),) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 15, hjust = 1, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  # symmetric axis per facet
  geom_blank(data = df_range, aes(x = "", y = max_val)) +
  geom_blank(data = df_range, aes(x = "", y = -max_val))
#print(p)

# --- save plot ---
ggsave(paste0("cams2_35_bis/plot/massdia_",expname1,"VS",expname2,"_",gsub("-","",sDate),"-",gsub("-","",eDate),".png"), 
       plot = p,
       width = 16*1.5, 
       height = 9*1.5, 
       units = "in",
       dpi = 300)

