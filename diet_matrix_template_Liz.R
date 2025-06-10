# build a diet matrix template for Liz
# groups
grps <- read.csv("data/GOA_Groups.csv", header = T)
pred_names <- grps %>% 
  filter(GroupType %in% c("MAMMAL","BIRD")) %>% 
  pull(Name)
pred_codes <- grps %>% 
  filter(GroupType %in% c("MAMMAL","BIRD")) %>% 
  pull(Code)

# amat
age_mat <- read.csv("data/age_at_mat_OY.csv", header = T)
age_mat <- age_mat %>% filter(Code %in% pred_codes)

# diet check (realized diets)
diet <- read.table("data/outputGOA01689_testDietCheck.txt", sep = " ", header  = T)

# handle diets
diet_liz <- diet %>%
  mutate(Time = Time / 365) %>%
  filter(Time > 45) %>% # focus on last 5 years of the run
  group_by(Predator, Cohort) %>%
  summarise(across(KWT:DR, mean)) %>% # mean (this is just to collapse the df)
  ungroup() %>%
  filter(Predator %in% pred_codes) %>%
  left_join(age_mat, by = c("Predator" = "Code")) %>% # bring in age at maturity
  mutate(stage = ifelse(Cohort < age_class_mat, "juvenile", "adult")) %>% #determine life stage
  group_by(Predator, stage, age_class_mat) %>%
  summarise(across(KWT:DR, mean)) %>% # mean (this is just to collapse the df)
  ungroup() 

# handle names
# will take pivoting twice
diet_liz <- diet_liz %>%
  pivot_longer(-c(Predator, stage, age_class_mat), names_to = 'Prey', values_to = 'Prop') %>%
  left_join((grps %>% select(Code, Name, LongName)), by = c('Predator'='Code')) %>%
  rename(Predator_Name = Name, Predator_LongName = LongName) %>%
  select(-Predator) %>%
  left_join((grps %>% select(Code, Name, LongName)), by = c('Prey'='Code')) %>%
  rename(Prey_Name = Name, Prey_LongName = LongName) %>%
  select(Predator_LongName, stage, age_class_mat, Prey_LongName, Prop) %>%
  pivot_wider(names_from = Prey_LongName, values_from = Prop)

# last fine details, such as more proper age at maturity
diet_liz <- diet_liz %>%
  left_join(grps %>% select(LongName, NumAgeClassSize),
            by = c("Predator_LongName" = "LongName")) %>%
  mutate(age_at_maturity = (age_class_mat) * (NumAgeClassSize) + 1) %>%
  select(Predator_LongName, stage, age_at_maturity, `Transient killer whales`:`Refractory detritus`) %>%
  arrange(Predator_LongName, desc(stage))

write.csv(diet_liz, "diet_matrix_template.csv", row.names = F)
