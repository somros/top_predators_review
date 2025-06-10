# Diets
# pull diets at equilibrium and make a series of plots
# turn it all into a quarto doc
library(tidyverse)
library(tidync)
library(ncdf4)

# load data
# groups
grps <- read.csv("data/GOA_Groups.csv", header = T)
pred_names <- grps %>% 
  filter(GroupType %in% c("MAMMAL","BIRD")) %>% 
  pull(Name)

# diet check (realized diets)
diet <- read.table("data/outputGOA01689_testDietCheck.txt", sep = " ", header  = T)
# nc file for WAA
nc <- nc_open("data/GOA_cb_summer.nc")
varnames <- names(nc$var)

# handle diets
diet_long <- diet %>%
  mutate(Time = Time / 365) %>%
  filter(Time > 45) %>% # focus on last 5 years of the run
  group_by(Predator, Cohort) %>%
  summarise(across(KWT:DR, mean)) %>%
  ungroup() %>%
  pivot_longer(-c(Predator, Cohort), names_to = 'Prey', values_to = 'Prop') %>%
  left_join((grps %>% select(Code, Name, LongName)), by = c('Predator'='Code')) %>%
  rename(Predator_Name = Name, Predator_LongName = LongName) %>%
  select(-Predator) %>%
  left_join((grps %>% select(Code, Name, LongName)), by = c('Prey'='Code')) %>%
  rename(Prey_Name = Name, Prey_LongName = LongName) %>%
  select(Prop, Predator_Name, Predator_LongName, Cohort, Prey_Name, Prey_LongName)%>%
  filter(Prop > 0.005)

diet_long <- diet_long %>% filter(Predator_Name %in% pred_names)

# Function to get fill value for a variable
get_fill_value <- function(nc, var_name) {
  fill_value <- ncatt_get(nc, var_name, "_FillValue")
  
  if (fill_value$hasatt) {
    return(fill_value$value)
  } else {
    return(NA)
  }
}


this_pred <- pred_names[1]

# make plots
p_diet <- diet_long %>%
  filter(Predator_Name == this_pred) %>%
  ggplot(aes(x = Cohort+1, y = Prop * 100, fill = Prey_LongName))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_x_continuous(breaks = 1:10)+
  #scale_fill_manual(values = colors)+
  theme_bw()+
  labs(x = '', y = "Diet preference (%)", fill = "Prey")
p_diet

# identify variables of interest
resn_vars <- varnames[grepl(this_pred, varnames) & grepl("ResN", varnames)]
structn_vars <- varnames[grepl(this_pred, varnames) & grepl("StructN", varnames)]

vars_df <- data.frame("var" = c(resn_vars, structn_vars)) %>%
  mutate(
    cohort = as.numeric(str_extract(var, "\\d+")),
    var_short = str_extract(var, "(ResN|StructN)$"),
    Name = str_remove(var, paste0(cohort, "_", var_short))
  )

vars_df$fillvalue <- sapply(vars_df$var, function(x) get_fill_value(nc, x))

# transform to kg
vars_df <- vars_df %>%
  mutate(weight_kg = fillvalue * 5.7 * 20 / 1000000)

# summarize
waa <- vars_df %>%
  group_by(Name, cohort) %>%
  summarise(waa = sum(weight_kg))

# close the nc file
nc_close(nc)
