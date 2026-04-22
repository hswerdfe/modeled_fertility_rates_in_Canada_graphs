library(scales)
library(cansim)
library(tidyverse)
library(ggplot2)
library(janitor)
library(lubridate)
library(glue)
library(feasts)
library(fpp3)
library(ggrepel)
library(magrittr)
library(viridis)
library(gganimate)
#install.packages("demography")
#library(demography)

library(tidyverse)
library(tsibble)
library(fable)






tlb_nm_fer <- '13-10-0418'
fer_dat_raw <- get_cansim(tlb_nm_fer)  |> clean_names()



theme_set(theme_minimal() +
            theme(
              axis.title = element_text(size = 15),
              panel.grid  = element_blank(),
              plot.title = element_text(hjust = 0.5, size = 30, color = 'darkgrey'),
              plot.subtitle = element_text(hjust = 0.5, size = 15, color = 'darkgrey')
            ))


dat_fert_all <-
  fer_dat_raw |> 
  filter(geo == 'Canada, place of residence of mother' &
           (
             str_detect(characteristics, 'Age-specific fertility rate, females') # |
             #str_detect(characteristics, 'Crude birth rate, live births per 1,000 population') 
           )
  ) |>
  mutate(yr = as.integer(ref_date)) |>
  mutate(age_group = str_squish(str_replace(str_remove(str_remove(characteristics, 'Age-specific fertility rate, females'), 'years'), '\\s+to\\s+','-'))) |>
  #filter(yr %in%   range(yr)) |>
  select(yr, value,   age_group) 



dat_fert_cohort <- 
  dat_fert_all |>
  mutate(age_group_min = as.integer(str_extract(age_group, '^([0-9]+)\\-([0-9]+)$',  group  = 1) ) ,
         age_group_max = as.integer(str_extract(age_group, '^([0-9]+)\\-([0-9]+)$',  group  = 2) )
  ) |>
  mutate(
    cohort_youngest = yr - age_group_min ,
    cohort_oldest = yr - age_group_max 
  ) |>
  mutate(
    cohort = glue('{cohort_oldest}-{cohort_youngest}')
  ) |>
  mutate(age_group = factor(fct_relevel(age_group, sort(unique((age_group)))), ordered = TRUE)) |>
  mutate(cohort = factor(fct_relevel(cohort, sort(unique((cohort)))), ordered = TRUE)) 



########################
# Use this for ALL Cohorts
#
all_cohorts <-
  tibble(cohort_youngest = seq(max(dat_fert_cohort$cohort_youngest), 1900, by = -5),) |>
  mutate(cohort_oldest = cohort_youngest - 4,
         cohort = glue('{cohort_oldest}-{cohort_youngest}')) |>
  distinct(cohort)


########################
# Use this for Cohorts, that are in most recent data set only
cohorts <- 
  dat_fert_cohort |>
  filter(yr == max(yr) ) |>
  distinct(cohort)

cohorts <- all_cohorts



p_dat <- 
  dat_fert_cohort |>
  filter(cohort %in% cohorts$cohort) 

# Get the maximum number of age groups available in the data
max_age_groups <- length(unique(p_dat$age_group))
# Find the youngest age group (e.g. "15-19")
youngest_age_group <- min(p_dat$age_group)


cum_fr_by_cohort <-
  dat_fert_cohort |>
  #filter(cohort %in% cohorts$cohort)  |>
  summarise(
    curr_fr = sum(value, na.rm = TRUE) * 5 / 1000,
    value = value,
    n_age_groups = n(),
    age_group = age_group  ,
    .by = c(cohort, yr)
  ) |>
  mutate(
    has_youngest = any(age_group == youngest_age_group),
    .by = c(cohort)
  )  |>
  # Drop cohorts that never observed the youngest age group
  #filter(has_youngest) |>
  group_by(cohort) |>
  arrange(yr, .by_group = TRUE) |>
  mutate(cum_fr = cumsum(curr_fr))


 
# Prepare the data as a tsibble
# We use 'yr' as the time index and 'age_group' as the key
ts_data <- cum_fr_by_cohort %>%
  ungroup() |>
  as_tsibble(index = yr, key = age_group)

# Fit a model to each age group
# Exponential Smoothing (ETS) or ARIMA are standard choices here
fit <- ts_data %>%
  model(
    # Log transformation ensures results stay > 0
    ets = ETS(log(value) ~ error("A") + trend("Ad") + season("N")),
    arima = ARIMA(log(value))
  )
report(fit) 
best_fit <- 
  accuracy(fit) |>
  filter(.model  == 'arima', .by = age_group)
  #filter(RMSE == min(RMSE), .by = age_group)

fc <- fit %>%
  forecast(h = "35 years")


predicted_measured_births <- 
  bind_rows(
    cum_fr_by_cohort |>
      ungroup() |> 
      select(age_group, yr, value) |>
      mutate(measured_predicted = 'measured'),
    fc |>
      as_tibble() |>
      filter(.model  == 'arima', .by = age_group) |>
      select (age_group, yr, .mean) |> mutate(measured_predicted = 'predicted') |>
      rename(value := .mean)
  ) 

pmb <- predicted_measured_births





# cum_fr_by_cohort |>
#   ungroup() |>
#   slice_max(order_by = value, by = yr ) |>
#   select(age_group, yr, value) |> View()


ggplot(mapping = aes(
  y        = value,
  x        = yr,
  group    = interaction(age_group, measured_predicted),
  color    = age_group,
  linetype = measured_predicted
)) +
geom_vline(xintercept = max(cum_fr_by_cohort$yr)+0.5) +
  # annotate(
  #   "text",
  #   x     = max(cum_fr_by_cohort$yr) + 1,
  #   y     = max(pmb$value)/2,          
  #   label = max(cum_fr_by_cohort$yr),
  #   angle = 90,
  #   #hjust = -0.1,        
  #   #vjust = -0.5,          
  #   size  = 5,
  #   color = "grey50",
  #   alpha = 0.5
  #)+
geom_point(data = pmb |> filter(measured_predicted == 'measured')) +
geom_line(data = pmb |> filter(measured_predicted == 'predicted')) +
geom_label_repel(
  #data = pmb |>filter(yr %in% range(yr) | yr == max(cum_fr_by_cohort$yr), .by = c(age_group)),
  data = pmb |>filter(yr %in% range(yr) , .by = c(age_group)),
  mapping = aes(label = glue("{round(value, 1)}")),
  #nudge_x = 1,
  #nudge_y = 5,
  alpha = 0.75,
  hjust = 0,
  
  #color = 'black',
  #fill = 'white',
  size =5,
  show.legend = FALSE
) +
  annotate(
    "text",
    x     = 2006.5,
    y     = 116,
    label = "In 2005, 30-34 Became The Most Fertile Age.",
    hjust = 0,
    size  = 4,
    color = "grey40"
  ) +
  annotate(
    "segment",
    x = 2006, xend = 2005,
    y = 114, yend = 101.0,
    arrow = arrow(length = unit(0.2, "cm")),
    color = "grey40"
  ) +
scale_color_viridis_d(option = 'H') +
#scale_color_viridis_d(option = 'B') +
  scale_x_continuous(
    breaks = seq(min(pmb$yr), max(pmb$yr), by = 5)
  ) +
  labs(
    x = glue(''), 
    y = glue('Births Per 1000 Females'),
    title = 'Birth Rates in Canada by Age Group.',
    subtitle = glue('Data After {max(cum_fr_by_cohort$yr)} is Predicted'),
    color = glue('Age Group of Mother'),
    caption = glue('Statscan Datatable: {tlb_nm_fer}')
  ) + 
  guides(
    color    = guide_legend(override.aes = list(shape = 16, linetype = 0, size = 3)),
    linetype = 'none',
    fill     = 'none'
  )



tfr_by_year <-
  pmb |>
  summarise(
    tfr = sum(value, na.rm = TRUE) * 5 / 1000,
    measured_predicted = 
      if (mean(measured_predicted == "measured", na.rm = TRUE) == 1){"measured"}else{'predicted'},
    .by = yr
  ) |>
  select(yr, tfr, measured_predicted) |>
  arrange(yr)

###################################
# Official Statscan TFR data, for comparison with our own calculations and predictions
# 
scenario_url = 'https://www150.statcan.gc.ca/n1/pub/17-20-0003/172000032026001-eng.htm'

##################
# Anchor points from StatCan Table 2
tfr_anchors <- tribble(
  ~year,  ~LG,  ~M1,  ~M2,  ~M3,  ~M4,  ~M5,  ~M6,  ~HG,  ~SA,  ~FA,
  2030,   1.12, 1.24, 1.24, 1.24, 1.24, 1.24, 1.24, 1.36, 1.36, 1.12,
  2050,   1.09, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.55, 1.55, 1.09,
  2075,   1.09, 1.32, 1.32, 1.32, 1.32, 1.32, 1.32, 1.55, 1.55, 1.09
)
scenarios <- c("LG", "M1", "M2", "M3", "M4", "M5", "M6", "HG", "SA", "FA")
scenarios_names <- c("Low Growth", "Medium Growth", "Medium Growth", "Medium Growth", "Medium Growth", "Medium Growth", "Medium Growth", "High Growth", "Slow Aging", "Fast Aging")


# Interpolate each scenario across all years 2025–2075
tfr_scenarios <- map(scenarios, \(scen) {
  anchors <- 
    tfr_anchors |> 
    select(year, tfr = all_of(scen)) |>
    bind_rows(
      tfr_by_year |> filter(measured_predicted == 'measured') |> slice_max(yr) |> select(yr, tfr) |> rename(year := yr)
    )
  
  tibble(year = 2025:2075) |>
    left_join(anchors, by = "year") |>
    mutate(tfr = approx(anchors$year, anchors$tfr, xout = year)$y,
           scenario = scen)
}) |>
  list_rbind() |>
  mutate(yr = year) |>
  select(yr, tfr, scenario) |>
  left_join(
    tibble(scenario = scenarios, scenario_names = scenarios_names),
    by = 'scenario'
  ) |>
  select(yr, tfr, scenario_names) |> distinct() |>
  filter(str_detect(scenario_names, 'Growth')) |>
  filter(yr <= tfr_by_year$yr  |> max()) |>
  mutate(measured_predicted = 'predicted')
  #filter(str_detect(scenario, 'M1|M2|M3|M4|M5|M6')) #|>
  #filter(str_detect(scenario, 'LG|HG')) 




# Improved Plotting Logic
tfr_by_year |>
  ggplot(aes(x = yr, y = tfr)) +
  # 1. Reference Line for Replacement Fertility
  #geom_hline(yintercept = 2.1, linetype = "dashed", color = "darkred", alpha = 0.3) +
  #annotate("text", x = mean(range(tfr_by_year$yr)), y = 2.15, label = "Replacement Level (2.1)", color = "darkred", size = 4) +
  
  # 2. Vertical Divider
  geom_vline(xintercept = 2024.5, color = "grey70") + 
  
  # 3. StatCan Scenarios (Background)
  geom_line(data = tfr_scenarios, aes(color = scenario_names), linewidth = 0.8, alpha = 0.6) +
  #geom_point(data = tfr_scenarios, aes(color = scenario_names), alpha = 0.6) +
  
  # 4. Measured Data (High Contrast)
  geom_line(data = filter(tfr_by_year, measured_predicted == 'measured'), linewidth = 1.2) +
  geom_point(data = filter(tfr_by_year, measured_predicted == 'measured'), size = 2) +
  
  # 5. Your ARIMA Model (Distinct Style)
  geom_line(data = filter(tfr_by_year, measured_predicted == 'predicted'), 
            color = "royalblue", linewidth = 1.2, linetype = "solid") +
  
  # 6. Final Touch: Labels
  geom_text_repel(
    data = tfr_scenarios |> filter(yr == max(yr)),
    aes(label = glue('Statscan\n{scenario_names} {round(tfr, 2)}'), color = scenario_names),
    nudge_x = 2, direction = "y", hjust = 0
  ) +
  # ARIMA specific model
  geom_label(
      data = tfr_by_year |> slice_max(yr),
      alpha=0.5,
      mapping = aes(
       x = yr + 1,
       y = tfr + 0.05,
       label = glue("ARIMA Model\nFertility, By Age Group {round(tfr, 2)}")
      ),
      fill = "royalblue", color = "white", fontface = "bold") +
  
  # Aesthetic tweaks
  scale_color_viridis_d(option = 'plasma', end = 0.8) +
  scale_x_continuous(
    breaks = seq(min(tfr_by_year$yr), max(tfr_by_year$yr), by = 5)
  ) +
  coord_cartesian(clip = 'off') + # Allows labels to bleed into margins
  #theme_minimal(base_size = 14) +
  theme(
    #plot.title = element_text(face = "bold", size = 22, color = "black"),
    plot.margin = margin(10, 50, 10, 10) # Room for labels on the right
  ) +
  labs(
    x = glue(''), 
    y = glue('Total Fertility Rate'),
    title = 'Total Fertility Rate in Canada Over Time',
    subtitle = glue('Data After {max(cum_fr_by_cohort$yr)} is Predicted'),
    caption = glue('Statscan Datatable: {tlb_nm_fer}, scenario : {scenario_url}'),
    color = glue('Official Stats')
  ) +
  guides(alpha = 'none', color = 'none', linetype = 'none')




cum_fr <- 
  pmb |> 
    mutate(age_group_min = as.integer(str_extract(age_group, '^([0-9]+)\\-([0-9]+)$',  group  = 1) ) ,
           age_group_max = as.integer(str_extract(age_group, '^([0-9]+)\\-([0-9]+)$',  group  = 2) )
    ) |>
    mutate(
      cohort_youngest = yr - age_group_min ,
      cohort_oldest = yr - age_group_max 
    ) |>
    mutate(
      cohort = glue('{cohort_oldest}-{cohort_youngest}')
    ) |>
    mutate(cohort = factor(fct_relevel(cohort, sort(unique((cohort)))), ordered = TRUE)) |>
    #filter(cohort %in% cohorts$cohort)  |>
  summarise(
    curr_fr = sum(value, na.rm = TRUE) * 5 / 1000,
    value = value,
    n_age_groups = n(),
    age_group = age_group  ,
    measured_predicted = 
      if (mean(measured_predicted == "measured", na.rm = TRUE) == 1){"measured"}else{'predicted'},   # 0–1 proportion
    .by = c(cohort, yr)
  ) |>
    mutate(
      has_youngest = any(age_group == youngest_age_group),
      .by = c(cohort)
    ) |>
    group_by(cohort) |>
    arrange(yr, .by_group = TRUE) |>
    mutate(cum_fr = cumsum(curr_fr))




cum_fr_partial <-
  cum_fr |> 
  filter(has_youngest) |>
  filter(cohort %in% cohorts$cohort) |>
  ungroup()

# Create bridge rows: last measured point for each cohort, relabelled as predicted
bridge <- cum_fr_partial |>
  filter(measured_predicted == "measured") |>
  slice_max(yr, by = cohort) |>
  mutate(measured_predicted = "predicted")

cum_fr_partial_connected <- bind_rows(cum_fr_partial, bridge) |>
  arrange(cohort, age_group)


cum_fr_partial_connected |>   
  ggplot(mapping = aes(
    x = as.integer(age_group), 
    y = cum_fr, 
    group = interaction(cohort, measured_predicted),
    color = as.character(cohort),
    linetype = measured_predicted,
    alpha = measured_predicted
  )) +
  geom_hline(yintercept = 2.1, color = "black") +
  geom_point(data = cum_fr_partial_connected |> filter(measured_predicted == 'measured')) +
  geom_line() +
  geom_label_repel(
    data = bind_rows(
      cum_fr_partial |> filter(age_group  %in% max(age_group ), .by = c(cohort))#,
      #cum_fr_partial_connected |> filter(measured_predicted  == 'measured') |> filter(yr  == max(yr), .by = c(cohort))
    ) |> distinct()

    ,
    mapping = aes(label = glue("{cohort}\n FR {round(cum_fr, 2)} {measured_predicted}")),
    nudge_x = 0.52,
    #nudge_y = 5,
    alpha = 1,
    hjust = 0.5,

    #color = 'black',
    #fill = 'white',
    size =3,
    show.legend = FALSE
  ) +
  #scale_color_viridis_d(option = 'B') +
  #scale_color_viridis_d(option = 'H', begin = 0.05, end = 0.85) +
  scale_color_manual(values = c(
    "1975-1979" = "#1b4f72",  # dark navy
    "1980-1984" = "#2980b9",  # blue
    "1985-1989" = "#16a085",  # teal
    "1990-1994" = "#e67e22",  # orange (replaces light green)
    "1995-1999" = "#c0392b",  # red
    "2000-2004" = "#8e44ad",  # purple
    "2005-2009" = "#2c3e50"   # near-black
  )) +
  scale_alpha_manual(values = c(measured = 1, predicted = 0.7)) +
  scale_x_continuous(breaks = sort(unique(as.integer(p_dat$age_group))), 
                     labels = \(.x){
                       levels(p_dat$age_group)[.x]
                     }) +
  labs(
    x = glue('Age of Mother when Giving Birth'), 
    y = glue('Cumulative Fertility Rate'),
    title = 'Completed Cohort Fertility in Canada',
    subtitle = glue('Data After {cum_fr |> filter(measured_predicted  == "measured") |> pull(yr) |> max()} (dashed line) is Predicted.'),
    color = glue('Birth/Cohort of Mother'),
    linetype = glue('Measured or Predicted'),
    caption = glue('Statscan Datatable: {tlb_nm_fer}')
  )+ 
  guides(#color = 'none',
         linetype = 'none', 
         fill = 'none',
         alpha = 'none'
  ) +
  theme(
    panel.grid.major.y = element_line(color = 'lightgrey',   linewidth  = 0.01, linetype  = 'solid'),
    panel.grid.minor.y = element_line(color = 'lightgrey',   linewidth  = 0.01, linetype  = 'solid'),
    panel.grid.major.x = element_line(color = 'lightgrey',   linewidth  = 0.01, linetype  = 'solid'),
    axis.text.x = element_text(angle = 0, size = 15,  hjust = 0.5, vjust = 0,  color = 'grey'),
    axis.text.y = element_text(size = 15, color = 'grey'),
    axis.title = element_text(size = 25, color = 'grey'),
    legend.title = element_text(angle = 0, size = 14,  hjust = 0.5, vjust = 0, ),
    plot.caption = element_text(color = 'grey', size = 10)
  ) +
  annotate(
    "text",
    x     = mean(range(as.integer(cum_fr_partial_connected$age_group))),
    y     = 2.0,
    label = 'Generational Replacement (2.1)',
    angle = 0,
    #hjust = -0.1,
    #vjust = -0.5,
    size  = 5,
    color = "black",
    alpha = 0.5
  )





p_dat_bg_lbl <- 
  p_dat |> mutate(
    age_group  = mean(range(as.integer(age_group))),
    value   = mean(range(value ))
  ) |>
  distinct(age_group, value, yr)

yr_rng <- range(p_dat$yr)

p_dat_cohort_lbl <- 
  p_dat |>
  mutate(lbl = glue("Mothers Birth Year\n{cohort}"))


anim <- 
  p_dat |>
  ggplot(aes(
    x = as.integer(age_group), 
    y = value
  )) +
  geom_text(
    data = p_dat_bg_lbl, 
    mapping = aes(
      label = as.character(yr), 
      x = as.integer(age_group), 
      y = value
    ), 
    inherit.aes = FALSE, 
    size = 50, 
    color = 'lightgrey',
    alpha = 0.5
  ) +
  geom_point(size = 4, aes(color = as.character(cohort))) + 
  geom_line(mapping = aes(
    color = as.character(cohort),
    group = cohort,
    linetype = as.character(cohort)
    
    
  ), 
  linewidth = 1.2)+ 
  # geom_label(
  #   data = cum_fr_by_cohort,
  #   mapping = aes(
  #     x = max(as.integer(p_dat$age_group)),   # pin to rightmost age group
  #     y = max(p_dat$value, na.rm = TRUE) * 0.9,  # near top of plot
  #     label = glue("TFR: {round(tfr, 2)}"),
  #     group = as.character(cohort),
  #     color = as.character(cohort)
  #   ),
  #   inherit.aes = FALSE,
  #   fill = 'white',
  #   alpha = 0.85,
  #   size = 5,
  #   hjust = 1
  # ) +
  geom_label(
    data = p_dat_cohort_lbl,
    mapping = aes(
      label = lbl,#glue("Mothers Birth Year\n{cohort}"),
      group =  as.character(cohort)
    ),
    nudge_x = 0.5,
    alpha = 0.75,
    hjust = 0.5,
    color = 'black',
    fill = 'white',
    size =5
  ) +
  scale_color_viridis_d(option = 'B') +
  #scale_color_viridis_d(option = 'H') +
  #scale_fill_viridis_d(option = 'H') +
  scale_fill_viridis_d(option = 'B') +
  scale_x_continuous(breaks = sort(unique(as.integer(p_dat$age_group))), 
                     labels = \(.x){
                       levels(p_dat$age_group)[.x]
                     }) +
  labs(
    x = glue('Age of Mother when Giving Birth'), 
    y = glue('Births Per 1000 Females'),
    title = 'Birth Rates in Canada by Cohort.',
    subtitle = glue('Data from {yr_rng[1]}-{yr_rng[2]}'),
    color = glue('Birth/Cohort of Mother'),
    linetype = glue('Birth/Cohort of Mother'),
    caption = glue('Statscan Datatable: {tlb_nm_fer}')
  ) + 
  guides(color = 'none',
         linetype = 'none', 
         fill = 'none'
  ) +
  theme(
    panel.grid.major.y = element_line(color = 'lightgrey',   linewidth  = 0.01, linetype  = 'solid'),
    panel.grid.minor.y = element_line(color = 'lightgrey',   linewidth  = 0.01, linetype  = 'solid'),
    panel.grid.major.x = element_line(color = 'lightgrey',   linewidth  = 0.01, linetype  = 'solid'),
    axis.text.x = element_text(angle = 0, size = 15,  hjust = 0.5, vjust = 0,  color = 'grey'),
    axis.text.y = element_text(size = 15, color = 'grey'),
    axis.title = element_text(size = 25, color = 'grey'),
    legend.title = element_text(angle = 0, size = 14,  hjust = 0.5, vjust = 0, ),
    plot.caption = element_text(color = 'grey', size = 10)
  ) +
  transition_reveal(yr)  +
  enter_fade() +
  exit_fade() +
  ease_aes('linear')


ap <- 
  animate(anim, 
          nframes = (length(unique(p_dat$yr))  * 40), 
          fps = 10, 
          end_pause = 40, 
          start_pause = 20,
          width = 1261,    # Set width in pixels
          height =  700,
          renderer = gifski_renderer()
  )
ap  
anim_save(file.path("cohort_birth_rates_by_age_and_year.gif"), 
          animation = ap)
