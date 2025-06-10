# produce new diets for top predators
library(tidyverse)
library(readxl)

grp <- read.csv("data/GOA_Groups.csv")
key <- grp %>% select(Code, LongName)

predcodes <- c("KWT", "WHT", "KWR", "DOL", "WHG", "WHH", "WHB", "SSL", "PIN", "BDF", "BDI", "BSF", "BSI")

# read in diets from Liz
newdiet <- read_xlsx("fromLiz/Atlantis mammal and seabird diets.xlsx",
                    sheet = 2)

# rename dolphins and ssl
newdiet$Predator_LongName <- gsub(" and porpoises", "", newdiet$Predator_LongName)
newdiet$Predator_LongName <- gsub("lions", "lion", newdiet$Predator_LongName)

# reshape
newdiet <- newdiet %>%
  select(-`...82`) %>%
  pivot_longer(-c(Species:Predator_LongName), names_to = "Prey_LongName", values_to = "Prop")

# aggregate by top predator weigthing by biomass contributions
newdiet <- newdiet %>%
  mutate(Prop = replace_na(Prop, 0)) %>%
  group_by(Predator_LongName, Prey_LongName) %>%
  summarize(Prop_mean = weighted.mean(Prop, `Proportion of biomass`)) %>%
  ungroup()

# now we need to turn this into a pprey matrix
pprey <- newdiet %>%
  left_join(key, by = c("Predator_LongName"="LongName")) %>%
  rename(PredCode = Code) %>%
  left_join(key, by = c("Prey_LongName"="LongName")) %>%
  rename(PreyCode = Code) %>%
  select(PredCode, PreyCode, Prop_mean) %>%
  mutate(PredCode = factor(PredCode, levels = predcodes),
         PreyCode = factor(PreyCode, levels = grp %>% pull(Code)))

pprey <- pprey[order(pprey$PredCode, pprey$PreyCode), ]


# reshape
pprey <- pprey %>%
  pivot_wider(names_from = PreyCode, values_from = Prop_mean)

# add sed
pprey <- pprey %>% mutate(DCSed = 0,
                          DLSed = 0,
                          DRSed = 0)

# write out pprey chunk
file_liz <- "liz_pprey.txt"
file.create(file_liz)
for(i in 1:length(predcodes)){
  for(j in 1:2){
    for(k in 1:2){
      lab <- paste0("pPREY",j,predcodes[i],k,"  81")
      vec <- pprey %>% filter(PredCode == predcodes[i]) %>% select(-PredCode) %>% paste(collapse = " ")
      cat(lab, file = file_liz, sep = "\n", append = T)
      cat(vec, file = file_liz, sep = "\n", append = T)
    }
  }
}
  