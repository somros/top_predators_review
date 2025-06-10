# calibrate predator diet based on difference between input and output
# will always need to comply with no values > 1

library(tidyverse)

#############################################
## Starting from 2091
## Here we use as target the PPREY matrix that contains Liz's values

target <- "liz_pprey.txt"
grps <- read.csv("data/GOA_Groups.csv")
codes <- c(grps$Code, c("DCsed","DLsed","DRsed"))

format_pprey <- function(file_path) {
  # Read the file
  lines <- readLines(file_path)
  
  # Remove empty lines
  lines <- lines[nchar(trimws(lines)) > 0]
  
  # Initialize result list
  result_list <- list()
  
  # Process pairs of lines (ID line + values line)
  for (i in seq(1, length(lines), by = 2)) {
    # First line contains PreyID and "81"
    id_line <- lines[i]
    id_elements <- unlist(strsplit(trimws(id_line), "\\s+"))
    prey_id <- id_elements[1]  # Extract just the PreyID part
    
    # Second line contains the 81 values
    if (i + 1 <= length(lines)) {
      values_line <- lines[i + 1]
      values <- as.numeric(unlist(strsplit(trimws(values_line), "\\s+")))
      
      # Combine PreyID with values
      result_list[[length(result_list) + 1]] <- c(prey_id, values)
    }
  }
  
  # Convert to data frame
  df <- as.data.frame(do.call(rbind, result_list), stringsAsFactors = FALSE)
  
  # Fix column types (first column as character, rest as numeric)
  df[,1] <- as.character(df[,1])
  df[,2:ncol(df)] <- lapply(df[,2:ncol(df)], as.numeric)
  
  # Set column names
  colnames(df) <- c("PreyID", codes)
  
  # handle the predator and stage identifiers
  df <- df %>%
    mutate(pred = str_extract(df[,1], "(?<=pPREY\\d)[A-Z]+(?=\\d$)"),
           prey_stage = as.numeric(gsub("pPREY","",substr(PreyID,1,6))),
           pred_stage = as.numeric(substr(PreyID, nchar(PreyID), nchar(PreyID)))) %>%
    select(pred, pred_stage, prey_stage, KWT:DRsed)
  
  # resgape
  df_long <- df %>%
    pivot_longer(-c(pred,pred_stage,prey_stage), names_to = "prey", values_to = "prop_target")
  
  return(df_long)
}

target_consumption <- format_pprey(target)

# now bring in the Atlantis output which has the realized consumption
# note that this file does not include prey cohort or life stage
run <- 2093
realized_consumption <- read.delim(paste0("data/outputGOA0",run,"_testDietCheck.txt"), sep = " ")
# handle time
ts <- unique(realized_consumption$Time)/365
max_t <- ts[(length(ts)-1)]*365
age_mat <- read.csv("data/age_at_mat_OY.csv")

realized_consumption <- realized_consumption %>%
  filter(Time == max_t) %>% # use the last time step but careful with truncated years that do not get written out
  filter(Predator %in% unique(target_consumption$pred)) %>% # focus on top preds here
  mutate(DCsed = 0, DLsed = 0, DRsed = 0) %>% # add sed columns
  left_join(age_mat, by = c("Predator"="Code")) %>% # go from cohort to life stage
  rowwise() %>%
  mutate(pred_stage = ifelse(Cohort<age_class_mat, 1, 2)) %>%
  ungroup() %>%
  select(Predator,pred_stage,KWT:DRsed) %>%
  pivot_longer(-c(Predator,pred_stage), names_to = "prey", values_to = "prop") %>%
  rename(pred = Predator,
         prop_realized = prop) %>%
  group_by(pred, pred_stage, prey) %>% # average props across age classes (now life stages)
  summarise(prop_realized = mean(prop_realized)) %>%
  ungroup()

# check props
realized_consumption %>% group_by(pred,pred_stage) %>% summarise(check = sum(prop_realized)) %>% pull(check) %>% summary()

# join
comp <- target_consumption %>%
  left_join(realized_consumption) %>%
  mutate(prop_ratio = prop_realized / prop_target)

# plot - how off are we?
comp %>%
  ggplot(aes(x = pred, y = prey, fill = prop_ratio))+
  geom_tile()+
  viridis::scale_fill_viridis()+
  facet_grid(~pred_stage)

# > summary(comp$prop_ratio)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   0.000   0.015   0.107   0.705   0.707  16.416    3404 

# tune - how to handle NA's?
new_pprey <- comp %>%
  rowwise() %>%
  mutate(new_pprey = ifelse(prop_target == 0, 0, prop_target / prop_ratio)) %>%
  ungroup() %>%
  mutate(new_pprey = ifelse(new_pprey>1,0.99,new_pprey)) %>%
  select(pred, pred_stage, prey, prey_stage, new_pprey)

# pivot wider
new_pprey_wide <- new_pprey %>%
  pivot_wider(names_from = prey, values_from = new_pprey)

# write out new pprey
predcodes <- unique(target_consumption$pred)
file_calib <- paste0("calibrated_pprey/pprey_calib_",run,".txt")
file.create(file_calib)
for(i in 1:length(predcodes)){
  for(j in 1:2){
    for(k in 1:2){
      lab <- paste0("pPREY",j,predcodes[i],k,"  81")
      vec <- new_pprey_wide %>% filter(pred == predcodes[i], prey_stage == j, pred_stage == k) %>% select(KWT:DRsed) %>% paste(collapse = " ")
      cat(lab, file = file_calib, sep = "\n", append = T)
      cat(vec, file = file_calib, sep = "\n", append = T)
    }
  }
}




# identify and index the PPREY matrix from the PRM file
diets_start <- grep("pPREY1KWT1", bio_prm) # flag of the first line in the PRM - change to first species
pprey_ind <- which(startsWith(x=bio_prm, "pPREY")==T)
diets_end <- max(pprey_ind)+2
DM_prm <- bio_prm[seq(diets_start,diets_end)]
names_pprey <- bio_prm[pprey_ind]
val_pprey <- bio_prm[pprey_ind+1]

FG <- gsub(" ","",unique(gsub("pPREY","",gsub('[[:digit:]]+', '', gsub("\\   .*","",names_pprey))))) # consumers
names_pprey_age <- unique(gsub("pPREY","",gsub('[[:digit:]]+', '',  gsub("\\   .*","", grep("pPREY1", names_pprey, value=T))))) # age-structured consumers

DM_to_format <- t(
  sapply(seq(1, length(val_pprey)),
         function(x){
           vec <- unlist(strsplit(val_pprey[x], " "))
           return(vec)
         })
)

for_order_FG <- FG_tab_Atl[,1] # groups from group file
colnames(DM_to_format) <- c(for_order_FG,c("DCsed","DLsed","DRsed"))

formatted_DM_original <- DM_to_format %>%
  as_tibble() %>%
  cbind(label=gsub("\\ .*","",names_pprey)) %>%
  cbind(PredatorCODE=gsub("pPREY","",gsub('[[:digit:]]+', '', gsub("\\ .*","",names_pprey)))) %>%
  cbind(PreyAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\ .*","", names_pprey)),1,1)%in%c(1,2),
                            substr(gsub("pPREY", "", gsub("\\ .*","", names_pprey)),1,1),
                            "1")) %>%  
  #    cbind(PreyAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\   .*","",names_pprey)),1,1)%in%c(1,2),
  #                                 substr(gsub("pPREY", "", gsub("\\   .*","",names_pprey)),1,1),
  #                                 "2")) %>%  
  cbind(PredatorAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\ .*","",names_pprey)),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))))%in%c(1,2), 
                                substr(gsub("pPREY", "", gsub("\\ .*","",names_pprey)),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey)))),
                                "2")
  ) %>% 
  mutate(PredatorAgeClass=ifelse(PredatorAgeClass==1,"Juvenile", "Adult"),
         PreyAgeClass=ifelse(PreyAgeClass==2,"Adult", "Juvenile")) %>%
  dplyr::select(c("label","PredatorCODE", "PreyAgeClass", "PredatorAgeClass", colnames(DM_to_format)))

# tis produces a long table with the initial PPREY values
# These are not necessarily proportional diets
formatted_DM_original_long <- formatted_DM_original %>%
  pivot_longer(5:ncol(formatted_DM_original), values_to="diet_init", names_to="Prey") %>%
  group_by(PredatorCODE, PredatorAgeClass, Prey) %>%
  mutate(diet_init= as.numeric(diet_init)) %>%
  summarize(diet_init=mean(diet_init, na.rm=T))%>%
  ungroup()

# These should add up to 0.1 (or close) for each "row"
formatted_DM_original_long %>%
  group_by(PredatorCODE, PredatorAgeClass) %>%
  summarize(check = sum(diet_init)) %>%
  pull(check) %>%
  hist()

#################################################################
#### Calibrate top predators and groundfish in the GOA model
#################################################################

### Read outputs
## Starting from version with amplified consumption, run 1689

# model output (the latest that we want to fix)
pred_diets <- read.table("C:/Users/Alberto Rovellini/Documents/GOA/Parametrization/output_files/data/out_1689/outputGOA01689_testDietCheck.txt", header=T)
# read in age at maturity - mind that this is not likely to align with the age at which the ontogenetic shift in diets occurs
age_mat <- read.table("../data/age_mat.txt", sep = ' ') %>%
  select(V1, V4) %>%
  purrr::set_names(c('Code','age')) %>%
  mutate(Code = gsub('_age_mat','',Code))

formatted_DM_pred_long <- pred_diets %>%
  filter(Time==365) %>% # can consider using later time steps
  pivot_longer(6:ncol(pred_diets), values_to="diet_prop", names_to="Prey") %>%
  left_join(age_mat, by = c('Predator'='Code')) %>%
  mutate(PredatorAgeClass=ifelse(!Predator %in% no_age,
                                 ifelse(Cohort < age, "Juvenile", "Adult"), "Adult")) %>% 
  group_by(Time, Predator, Prey, PredatorAgeClass) %>% 
  summarize(diet_prop_est=sum(diet_prop, na.rm=T)) %>% # summing contributions to cohorts into life stages
  group_by(Time, Predator, PredatorAgeClass) %>%
  mutate(diet_prop_tot=sum(diet_prop_est, na.rm=T)) %>% 
  mutate(diet_prop_est=diet_prop_est/diet_prop_tot) %>% # rescales the contributions by the total, to make sure there are proportional (do we want this though?)
  ungroup() %>%
  dplyr::select(-diet_prop_tot, -Time) #%>%
#drop_na() # CAREFUL: dropping all biomass pools

# # write out matrices for Claude
# write.table(formatted_DM_original_long %>%
#               ungroup() %>%
#               select(PredatorCODE,PredatorAgeClass,Prey,diet_init)%>%
#               mutate(diet_init = diet_init*10),#%>%
#               #filter(PredatorCODE=="ATF"), 
#             "matrixA.csv", 
#             row.names = F,
#             col.names = c("Predator","PredStage","Prey","Value"),
#             sep=",")
# write.table(formatted_DM_pred_long %>%
#               select(Predator,PredatorAgeClass,Prey,diet_prop_est),#%>%
#               #filter(Predator=="ATF"), 
#             "matrixC.csv", 
#             row.names = F,
#             col.names = c("Predator","PredStage","Prey","Value"),
#             sep=",")

# IMPORTANT: this last step above rescales the PPREY contributions so that they sum up to 1.
# In the GOA model we needed to increase individual PPREY values, and the total did not sum up to 1
# Rescaling PPREY entries this way may help with proportional diet contributions, but if may also deflate consumption
# The PPREY matrix is defined as an availability matrix, and there is no real requirement that the values entered are proportional
# We may want to skip this step, let's try with and without rescaling

# check props here:
tt <- formatted_DM_pred_long %>%
  group_by(Predator, PredatorAgeClass) %>%
  summarise(check  = sum(diet_prop_est))

# calculate correction factors
corr_diet <- formatted_DM_original_long %>% # takes original diets (i.e., true diets we aim for)
  group_by(label, PredatorCODE, PredatorAgeClass) %>% # keeping the label is what keeps the prey stage
  mutate(sum_diet=sum(diet_init, na.rm=T)) %>% # sum of original diets, adds up to close to 0.1 for us
  ungroup() %>%
  mutate(pc_diet_init=diet_init/sum_diet) %>% # makes proportions (0-1) from original PPREY values
  left_join(formatted_DM_pred_long, by=c("PredatorCODE"="Predator", "Prey", "PredatorAgeClass")) %>% # joins with output consumptions
  group_by(label, PredatorCODE, PredatorAgeClass) %>%
  mutate(tot_diet_prop_est=sum(diet_prop_est, na.rm=T)) %>% # this seems redundant with the reapportioning above
  ungroup() %>%
  mutate(diet_prop_est=diet_prop_est/tot_diet_prop_est) %>%
  dplyr::select(-tot_diet_prop_est) %>%
  mutate(ratio_diet=diet_prop_est/pc_diet_init) %>% # this divides estimated/expected, aka Atlantis/true. It works on PROPORTIONS though
  mutate(ratio_diet=ifelse(diet_init==0, 0, ratio_diet)) %>%
  mutate(ratio_diet=ifelse(diet_prop_est==0, 1, ratio_diet))

# this step returned a set of multipliers for each entry of the PPREY matrix
# Now these must be applied to the most recent PPREY matrix to make it match the intended diets

# now package it for the groups of interest.
# some examples:
top_preds <- FG_tab_Atl %>%
  filter(GroupType %in% c("MAMMAL","BIRD","SHARK")) %>%
  pull(Code)
tier3 <- c("ATF","COD","POL","ATF","POP","SBF","FFS","FFD","RFS","HAL","REX","FHS")
forage <- c("HER","EUL","CAP","SAN","FOL")

todo <- tier3

prio_corr <- corr_diet %>%
  filter(PredatorCODE %in% todo) %>%
  filter(ratio_diet!=0) %>% # we do not want 0's?
  #filter(ratio_diet!=1) %>% # this is like doing nothing (i.e., original and predicted match)
  group_by(label, PredatorCODE, PredatorAgeClass, Prey) %>%
  summarize(mean_corr=mean(ratio_diet, na.rm=T))#, # this step seems redundant for GOA (hence why SD = NA, the grouping results in the sma enumber of rows)
#var_corr=sd(ratio_diet, na.rm=T))

# prio_corr %>%
#   ggplot() + 
#   geom_bar(aes(x=PredatorCODE, y=1/mean_corr, fill=PredatorAgeClass), position="dodge", stat="identity") +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=45)) +
#   facet_wrap(~Prey, scales="free_y")

# These need to be applied to the "recent" PPREY matrix as PPREY/CORR, because they were computed as predicted/true
# flip CORR so that it becomes multiplicative
# NOTE: # some of these are really high
# What this means is that the resulting PPREY will be really high, maybe close to 1
# This may be warranted in some cases but not in others, and it is difficult to generalize
# For now, ignore correction on Pteropods (we will have to fix this group but for the time being leave them)
prio_corr <- prio_corr %>%
  mutate(corr = 1/mean_corr) %>%
  select(-mean_corr) %>%
  mutate(corr = ifelse(Prey == "PTE", 1, corr))

# This bit seems conceptually different from what I need to do
# It is probably because their "true" and "latest" diet matrix were the same
# This is not the case for the GOA - the two are now really different
# I want to correct the current diet matrix, not the original
# Comment out but leave it here for future reference for now
###################################################################################
# change_diet <- corr_diet %>%
#   filter(PredatorCODE %in% todo) %>%
#   mutate(ratio_diet=ifelse(is.na(ratio_diet), 1, ratio_diet))
# 
# calc_new_pprey  <- formatted_DM_original_long %>% # take true diets
#   left_join(change_diet[, c("label","PredatorCODE","PredatorAgeClass", "Prey", "ratio_diet")], # join with correction
#             by=c("label","PredatorCODE","PredatorAgeClass", "Prey")) %>%
#   group_by(label) %>%
#   mutate(sum_orig=sum(diet_init, na.rm=T)) %>% # get total initial sum of PPREY
#   ungroup() %>%
#   mutate(new_diet=diet_init/ratio_diet) %>% # this takes the true diet and divides it by the correction factor?
#   group_by(label) %>%
#   mutate(sum_new=sum(new_diet, na.rm=T)) %>%
#   ungroup() %>%
#   mutate(sum_new=ifelse(sum_new==0, 1, sum_new)) %>%
#   mutate(verif=new_diet/sum_new) %>%
#   mutate(new_diet=ifelse(is.nan(new_diet), 0, new_diet)/sum_new*sum_orig)
# 
# readytransfo <- calc_new_pprey  %>%
#   dplyr::select(label, Prey, new_diet) %>%
#   pivot_wider(values_from=new_diet, names_from=Prey) %>%
#   dplyr::select(all_of(label,for_order_FG,DCsed,DLsed,DRsed))
###################################################################################
change_diet <- prio_corr

# read in new diet matrix (make this a function as it is the same workflow used above)
run <- 1689
latestAtlmat <- paste0("../../output_files/data/out_",run,"/GOAbioparam_test.prm") #PRM with recent diets to calibrate
bio_prm_latest <- readLines(latestAtlmat)

# identify and index the PPREY matrix from the PRM file
diets_start <- grep("pPREY1KWT1", bio_prm_latest) # flag of the first line in the PRM - change to first species
pprey_ind <- which(startsWith(x=bio_prm_latest, "pPREY")==T)
diets_end <- max(pprey_ind)+2
DM_prm <- bio_prm_latest[seq(diets_start,diets_end)]
names_pprey <- bio_prm_latest[pprey_ind]
val_pprey <- bio_prm_latest[pprey_ind+1]

FG <- gsub(" ","",unique(gsub("pPREY","",gsub('[[:digit:]]+', '', gsub("\\   .*","",names_pprey))))) # consumers
names_pprey_age <- unique(gsub("pPREY","",gsub('[[:digit:]]+', '',  gsub("\\   .*","", grep("pPREY1", names_pprey, value=T))))) # age-structured consumers

DM_to_format_latest <- t(
  sapply(seq(1, length(val_pprey)),
         function(x){
           vec <- unlist(strsplit(val_pprey[x], " "))
           return(vec)
         })
)

colnames(DM_to_format_latest) <- c(for_order_FG,c("DCsed","DLsed","DRsed"))

formatted_DM_latest <- DM_to_format_latest %>%
  as_tibble() %>%
  cbind(label=gsub("\\ .*","",names_pprey)) %>%
  cbind(PredatorCODE=gsub("pPREY","",gsub('[[:digit:]]+', '', gsub("\\ .*","",names_pprey)))) %>%
  cbind(PreyAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\ .*","", names_pprey)),1,1)%in%c(1,2),
                            substr(gsub("pPREY", "", gsub("\\ .*","", names_pprey)),1,1),
                            "1")) %>%  
  cbind(PredatorAgeClass=ifelse(substr(gsub("pPREY", "", gsub("\\ .*","",names_pprey)),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))))%in%c(1,2), 
                                substr(gsub("pPREY", "", gsub("\\ .*","",names_pprey)),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey))),
                                       nchar(gsub("pPREY", "", gsub("\\ .*","",names_pprey)))),
                                "2")
  ) %>% 
  mutate(PredatorAgeClass=ifelse(PredatorAgeClass==1,"Juvenile", "Adult"),
         PreyAgeClass=ifelse(PreyAgeClass==2,"Adult", "Juvenile")) %>%
  dplyr::select(c("label","PredatorCODE", "PreyAgeClass", "PredatorAgeClass", colnames(DM_to_format)))

# tis produces a long table with the initial PPREY values
# These are not necessarily proportional diets
formatted_DM_latest_long <- formatted_DM_latest %>%
  pivot_longer(5:ncol(formatted_DM_latest), values_to="diet_latest", names_to="Prey") %>%
  mutate(diet_latest= as.numeric(diet_latest)) %>%
  ungroup()

# for Claude
# write.table(formatted_DM_latest_long %>%
#               group_by(PredatorCODE, PredatorAgeClass, Prey) %>%
#               summarize(diet_latest=mean(diet_latest, na.rm=T)) %>% 
#               ungroup() %>%
#               select(PredatorCODE,PredatorAgeClass,Prey,diet_latest),#%>%
#               #filter(PredatorCODE=="ATF"), 
#             "matrixB.csv", 
#             row.names = F,
#             col.names = c("Predator","PredStage","Prey","Value"),
#             sep=",")

# now apply the correction factors
calc_new_pprey  <- formatted_DM_latest_long %>% # take true diets
  filter(PredatorCODE %in% todo) %>%
  left_join(change_diet[, c("label","PredatorCODE","PredatorAgeClass", "Prey", "corr")], # join with correction - this causes most of corr to be expanded as NA (they did not appear in the diets?)
            by=c("label","PredatorCODE","PredatorAgeClass", "Prey")) %>%
  group_by(label) %>%
  mutate(sum_latest=sum(diet_latest, na.rm=T)) %>% # get total sum of PPREY in the latest matrix
  ungroup() %>%
  mutate(new_diet=diet_latest * corr) %>% # this takes the latest diet and multiplies it by the correction factor
  group_by(label) %>%
  mutate(sum_new=sum(new_diet, na.rm=T)) %>% # some of these will be enormous
  ungroup() %>%
  mutate(sum_new=ifelse(sum_new==0, 1, sum_new)) %>%
  mutate(new_pprey=ifelse(is.nan(new_diet), 0, new_diet)/sum_new) # this reproportions the new_diet term and calculates the final pprey

#Check that the new diets add up to 1
calc_new_pprey %>%
  mutate(new_pprey = replace_na(new_pprey, 0)) %>%
  group_by(label) %>%
  summarize(check = sum(new_pprey)) %>%
  pull(check)

# view (not sure if helpful)
calc_new_pprey %>%
  mutate(new_pprey = replace_na(new_pprey, 0)) %>%
  mutate(PreyAgeClass = substr(label, 6, 6)) %>%
  mutate(PreyAgeClass = ifelse(PreyAgeClass == 1, "Juvenile", "Adult")) %>%
  ggplot(aes(x = Prey, y = PredatorCODE, fill = new_pprey))+
  geom_tile()+
  scale_fill_viridis()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(PreyAgeClass ~ PredatorAgeClass)

# turn into a table to be transformed into the pprey matrix
newpprey <- calc_new_pprey  %>%
  mutate(new_pprey = replace_na(new_pprey, 0)) %>%
  dplyr::select(label, Prey, new_pprey) %>%
  pivot_wider(values_from=new_pprey, names_from=Prey) %>%
  dplyr::select(label,all_of(for_order_FG),DCsed,DLsed,DRsed)

# pull the modified groups from here
pprey_file <- paste0('calibrated_',run,todo,'.txt')
file.create(pprey_file)

for(i in 1:length(names_pprey)){
  
  this_pprey_name <-  gsub('(81).*','\\1',names_pprey[i])
  this_label <- gsub('\\  81.*','',names_pprey[i])
  
  if(this_label %in% unique(newpprey$label)){
    
    cat(this_pprey_name, file=pprey_file, append=TRUE,'\n\n')
    cat(unlist(newpprey %>% filter(label == this_label) %>% select(-label)), file=pprey_file, append=TRUE, '\n')
    
  } else {
    
    cat(this_pprey_name, file=pprey_file, append=TRUE,'\n\n')
    cat(val_pprey[i], file=pprey_file, append=TRUE, '\n')
    
  }
  
}
