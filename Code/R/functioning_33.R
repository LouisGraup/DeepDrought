######################
###### SETTINGS ######
######################


# set start and end date

start_date <- as.POSIXct("2024-05-03 00:00:00")
end_date <- as.POSIXct("2024-09-29 23:59:00")






### creat folder for figures

dir.create(paste0("figures_", Sys.Date()))


# 0. libraries and data ---------------------------------------------------
# section about connection to the database

library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rLWFpg)
library(patchwork)


conn <- db_connect('grauplou', rstudioapi::askForPassword("Database password"), db_host = 'pgdblwf', db_name = 'lwf')

# 1. get data --------------------------------------------------

db_tbl(conn, 'ada', 'v_messvar', retrieve = FALSE)

metadata.df <- db_tbl(conn, 'ada', 'v_messvar', retrieve = TRUE) %>%
  dplyr::rename_with(toupper) %>%
  filter(PROJECT_NAME %in% "Pfynwald irrigation",
         (grepl("PFY_SC", INSTALLATION_NAME)),
         VARNAME_NAME %in% c("Atmospheric vapour pressure deficit", "Valve normally closed", "Water pressure", "Wind speed", "Wind direction", "Radiation shortwave incoming", "Photosynthetic photon flux density"),
         VARFREQ_NAME %in% '1 Minutes')

# retrieve the data from DB
meteo.df <- ada_get_data(conn,
                         installation = unique(metadata.df$INSTALLATION_ID),
                         messvar = unique(metadata.df$MESSVAR_ID),
                         messtime_from = start_date,
                         messtime_to = end_date,
                         #faults_correct = TRUE,
                         retrieve = TRUE) %>%
  dplyr::rename_with(toupper)

# Add additional columns from metadata.df to meteo.df by joining on MESSVAR_ID
meteo.df <- merge(meteo.df, metadata.df[, c("MESSVAR_ID", "MESSVAR_NAME", "INSTALLATION_NAME", "TREATMENT")], by = "MESSVAR_ID", all.x = TRUE)



# Calculate the median for the reference VPD of each MESSVAR_NAME and TREATMENT for each time step (MESSTIME)
median_values <- meteo.df %>%
  group_by(MESSTIME, MESSVAR_NAME, TREATMENT) %>%
  summarise(MEDIAN_REF = median(MESSVAL, na.rm = TRUE), .groups = 'drop')
# Join the median values back to meteo.df
meteo.df <- meteo.df %>%
  left_join(median_values, by = c("MESSTIME", "MESSVAR_NAME", "TREATMENT"))



# Step 1: Filter data for MESSVAR_NAME = "VPD_ic_ce" and TREATMENT in "irrigation_vpd" or "roof_vpd"
filtered_meteo.df <- meteo.df %>%
  filter(MESSVAR_NAME == "VPD_ic_ce", TREATMENT %in% c("irrigation_vpd", "roof_vpd"))
# Step 2: Calculate the median of MESSVAL for each MESSTIME
median_values <- filtered_meteo.df %>%
  group_by(MESSTIME) %>%
  summarise(MEDIAN_VPD = median(MESSVAL, na.rm = TRUE), .groups = 'drop')
# Step 3: Join the median values back to meteo.df
meteo.df <- meteo.df %>%
  left_join(median_values, by = c("MESSTIME"))


# fill MEDIAN_REF of reference (median of control plots) to entire time step
meteo.df <- meteo.df %>%
  group_by(MESSTIME) %>%
  mutate(
    MEDIAN_REF = unique(MEDIAN_REF[MESSVAR_NAME == "VPD_ic_ce" & TREATMENT == "control"])
  ) %>%
  ungroup()



# get absolute VPD change
meteo.df$DELTA_ABS <- meteo.df$MESSVAL-meteo.df$MEDIAN_REF
# # get absolute VPD reduction
# meteo.df$DELTA_ABS <- meteo.df$MEDIAN_REF-meteo.df$MESSVAL

# get relative VPD reduction
meteo.df$DELTA_REL <- meteo.df$DELTA_ABS/meteo.df$MEDIAN_REF*100
# # get relative VPD reduction
# meteo.df$DELTA_REL <- meteo.df$DELTA_ABS/meteo.df$MEDIAN_REF*100



### get time steps, when system was operating (based on water pressure and valve state)

###### define operation time of the system by water pressure
# filter for water pressure
filtered_data <- meteo.df %>%
  filter(MESSVAR_NAME == "WaterPressure_Avg")
# set time 
median_values <- filtered_data %>%
  group_by(MESSTIME) %>%
  summarise(MEDIAN_WATER = median(MESSVAL, na.rm = TRUE), .groups = 'drop')
# get time steps where median above threshold of 60bar
water_pressure <- median_values %>%
  filter(MEDIAN_WATER > 60) %>%
  pull(MESSTIME)
###### define operation time of the system by valve1 state
# filter for valve state
filtered_data <- meteo.df %>%
  filter(MESSVAR_NAME == "Valve_1_Avg")
# set time 
median_values <- filtered_data %>%
  group_by(MESSTIME) %>%
  summarise(MEDIAN_WATER = median(MESSVAL, na.rm = TRUE), .groups = 'drop')
# get time steps where median above threshold
valve1_state <- median_values %>%
  filter(MEDIAN_WATER < -0.99) %>%
  pull(MESSTIME)
###### define operation time of the system by valve2 state
# filter for valve state
filtered_data <- meteo.df %>%
  filter(MESSVAR_NAME == "Valve_2_Avg")
# set time 
median_values <- filtered_data %>%
  group_by(MESSTIME) %>%
  summarise(MEDIAN_WATER = median(MESSVAL, na.rm = TRUE), .groups = 'drop')
# get time steps where median above threshold
valve2_state <- median_values %>%
  filter(MEDIAN_WATER < -0.99) %>%
  pull(MESSTIME)
###### get time steps with water pressure above threshold, and valve1 and/or valve2 open
time_operating <- water_pressure[which(water_pressure %in% unique(c(valve1_state, valve2_state)))]

###### add operation status to meteo.df_complete
# Create the OPERATION_STATUS column, initialized with 0
meteo.df_complete <- meteo.df
meteo.df_complete$OPERATION_STATUS <- 0
# Update the OPERATION_STATUS to 1 where MESSTIME matches a value in time_operating
meteo.df_complete$OPERATION_STATUS[meteo.df_complete$MESSTIME %in% time_operating] <- 1

###### select time steps from meteo.df when system was operational
meteo.df <- meteo.df %>%
  filter(MESSTIME %in% time_operating)




# ###### define operation time of the system by difference between reference VPD and manipulated VPD in ic_ce
# # Step 1: Filter the data for MESSVAR_NAME == "VPD_ic_ce" and TREATMENT == "control"
# test <- meteo.df %>%
#   filter(MESSVAR_NAME == "VPD_ic_ce", TREATMENT == "control") %>%
#   mutate(DIFF = ifelse(MEDIAN_REF - MEDIAN_VPD >= 0.1, 999, 0)) #%>%
#   #mutate(DIFF = ifelse(MEDIAN_REF - MEDIAN_VPD >= MEDIAN_REF*0.2, 999, 0))
# # Step 2: Calculate the percentage of time with system operation (DIFF == 999)
# percentage_system_operation <- mean(test$DIFF == 999) * 100
# print(paste0("% time of system operation: ", percentage_system_operation))
# # Step 3: Filter meteo.df for rows where MESSTIME is in test$MESSTIME with DIFF == 999
# meteo.df <- meteo.df %>%
#   filter(MESSTIME %in% test$MESSTIME[test$DIFF == 999])


# ### get percentiles of observed reference VPD
# percentile_90 <- quantile(meteo.df$MEDIAN_REF, 0.90, na.rm = TRUE)
# meteo.df <- meteo.df[which(meteo.df$MEDIAN_REF >= percentile_90), ]




###### figures front to back


### figure front to back relative change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("VPD_ic", MESSVAR_NAME) & grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("fl", MESSVAR_NAME) ~ "front",
    grepl("fr", MESSVAR_NAME) ~ "front",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("bl", MESSVAR_NAME) ~ "back",
    grepl("br", MESSVAR_NAME) ~ "back",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("front", "center", "back"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_REL, 0.25)
Q3 <- quantile(filtered_data$DELTA_REL, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/FrontBack_rel_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_REL, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (%)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("front" = "blue", "center" = "deepskyblue", "back" = "cyan"))  # Custom color palette

dev.off()



### figure front to back absolute change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("VPD_ic", MESSVAR_NAME) & grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("fl", MESSVAR_NAME) ~ "front",
    grepl("fr", MESSVAR_NAME) ~ "front",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("bl", MESSVAR_NAME) ~ "back",
    grepl("br", MESSVAR_NAME) ~ "back",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("front", "center", "back"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_ABS, 0.25)
Q3 <- quantile(filtered_data$DELTA_ABS, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/FrontBack_abs_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_ABS, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (kPa)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("front" = "blue", "center" = "deepskyblue", "back" = "cyan"))  # Custom color palette

dev.off()






###### figures left to right


### figure front to back relative change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("VPD_ic", MESSVAR_NAME) & grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("fl", MESSVAR_NAME) ~ "left",
    grepl("bl", MESSVAR_NAME) ~ "left",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("fr", MESSVAR_NAME) ~ "right",
    grepl("br", MESSVAR_NAME) ~ "right",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("left", "center", "right"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_REL, 0.25)
Q3 <- quantile(filtered_data$DELTA_REL, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/LeftRight_rel_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_REL, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (%)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("left" = "red", "center" = "yellow", "right" = "deepskyblue"))  # Custom color palette

dev.off()


### figure front to back absolute change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("VPD_ic", MESSVAR_NAME) & grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("fl", MESSVAR_NAME) ~ "left",
    grepl("bl", MESSVAR_NAME) ~ "left",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("fr", MESSVAR_NAME) ~ "right",
    grepl("br", MESSVAR_NAME) ~ "right",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("left", "center", "right"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_ABS, 0.25)
Q3 <- quantile(filtered_data$DELTA_ABS, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/LeftRight_abs_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_ABS, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (kPa)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("left" = "red", "center" = "yellow", "right" = "deepskyblue"))  # Custom color palette

dev.off()



###### figures single positions


### figure front to back relative change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("VPD_ic", MESSVAR_NAME) & grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("fl", MESSVAR_NAME) ~ "front_left",
    grepl("fr", MESSVAR_NAME) ~ "front_right",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("bl", MESSVAR_NAME) ~ "back_left",
    grepl("br", MESSVAR_NAME) ~ "back_right",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("front_left", "front_right", "center", "back_left", "back_right"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_REL, 0.25)
Q3 <- quantile(filtered_data$DELTA_REL, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/AllPos_rel_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 4, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_REL, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (%)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("front_left" = "red", "front_right" = "tomato", "center" = "yellow", "back_left" = "dodgerblue2", "back_right" = "deepskyblue"))  # Custom color palette

dev.off()



### figure front to back absolute change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("VPD_ic", MESSVAR_NAME) & grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("fl", MESSVAR_NAME) ~ "front_left",
    grepl("fr", MESSVAR_NAME) ~ "front_right",
    grepl("ce", MESSVAR_NAME) ~ "center",
    grepl("bl", MESSVAR_NAME) ~ "back_left",
    grepl("br", MESSVAR_NAME) ~ "back_right",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("front_left", "front_right", "center", "back_left", "back_right"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_ABS, 0.25)
Q3 <- quantile(filtered_data$DELTA_ABS, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/AllPos_abs_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 4, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_ABS, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (kPa)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("front_left" = "red", "front_right" = "tomato", "center" = "yellow", "back_left" = "dodgerblue2", "back_right" = "deepskyblue"))  # Custom color palette

dev.off()




###### figures above to below canopy


### figure front to back relative change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("PFY_SC14", INSTALLATION_NAME),
         grepl("VPD", MESSVAR_NAME))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("oc", MESSVAR_NAME) ~ "above_canopy",
    grepl("uc", MESSVAR_NAME) ~ "below_canopy",
    grepl("ic_fl", MESSVAR_NAME) ~ "in_canopy",
    grepl("ic_fr", MESSVAR_NAME) ~ "in_canopy",
    grepl("ic_bl", MESSVAR_NAME) ~ "in_canopy",
    grepl("ic_br", MESSVAR_NAME) ~ "in_canopy",
    grepl("ce", MESSVAR_NAME) ~ "center",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("above_canopy", "in_canopy", "center", "below_canopy"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_REL, 0.25)
Q3 <- quantile(filtered_data$DELTA_REL, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/AboveBelow_rel_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3.5, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_REL, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (%)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("above_canopy" = "deepskyblue", "in_canopy" = "lightgreen", "center" = "yellow", "below_canopy" = "brown"))  # Custom color palette

dev.off()



### figure front to back absolute change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("PFY_SC14", INSTALLATION_NAME),
         grepl("VPD", MESSVAR_NAME))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("oc", MESSVAR_NAME) ~ "above_canopy",
    grepl("uc", MESSVAR_NAME) ~ "below_canopy",
    grepl("ic_fl", MESSVAR_NAME) ~ "in_canopy",
    grepl("ic_fr", MESSVAR_NAME) ~ "in_canopy",
    grepl("ic_bl", MESSVAR_NAME) ~ "in_canopy",
    grepl("ic_br", MESSVAR_NAME) ~ "in_canopy",
    grepl("ce", MESSVAR_NAME) ~ "center",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("above_canopy", "in_canopy", "center", "below_canopy"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_ABS, 0.25)
Q3 <- quantile(filtered_data$DELTA_ABS, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/AboveBelow_abs_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3.5, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_ABS, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (kPa)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("above_canopy" = "deepskyblue", "in_canopy" = "lightgreen", "center" = "yellow", "below_canopy" = "brown"))  # Custom color palette

dev.off()



###### figures above to below canopy individual measurement points


### figure front to back relative change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("PFY_SC14", INSTALLATION_NAME),
         grepl("VPD", MESSVAR_NAME))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("oc_fl", MESSVAR_NAME) ~ "above_canopy_front_left",
    grepl("oc_fr", MESSVAR_NAME) ~ "above_canopy_front_right",
    grepl("oc_bl", MESSVAR_NAME) ~ "above_canopy_back_left",
    grepl("oc_br", MESSVAR_NAME) ~ "above_canopy_back_right",
    grepl("uc_fl", MESSVAR_NAME) ~ "below_canopy_front_left",
    grepl("uc_fr", MESSVAR_NAME) ~ "below_canopy_front_right",
    grepl("uc_bl", MESSVAR_NAME) ~ "below_canopy_back_left",
    grepl("uc_br", MESSVAR_NAME) ~ "below_canopy_back_right",
    grepl("ic_fl", MESSVAR_NAME) ~ "in_canopy_front_left",
    grepl("ic_fr", MESSVAR_NAME) ~ "in_canopy_front_right",
    grepl("ic_bl", MESSVAR_NAME) ~ "in_canopy_back_left",
    grepl("ic_br", MESSVAR_NAME) ~ "in_canopy_back_right",
    grepl("ce", MESSVAR_NAME) ~ "center",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("above_canopy_front_left", "above_canopy_front_right",
                                                         "above_canopy_back_left", "above_canopy_back_right",
                                                         "in_canopy_front_left", "in_canopy_front_right",
                                                         "in_canopy_back_left", "in_canopy_back_right",
                                                         "center",
                                                         "below_canopy_front_left", "below_canopy_front_right",
                                                         "below_canopy_back_left", "below_canopy_back_right"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_REL, 0.25)
Q3 <- quantile(filtered_data$DELTA_REL, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/AboveBelowInd_rel_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 6, height = 4)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_REL, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (%)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) #+  # Dynamically set y-axis limits to exclude outliers
# scale_fill_manual(values = c("above_canopy" = "deepskyblue", "in_canopy" = "lightgreen", "center" = "yellow", "below_canopy" = "brown"))  # Custom color palette

dev.off()



### figure front to back absolute change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("PFY_SC14", INSTALLATION_NAME),
         grepl("VPD", MESSVAR_NAME))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("oc_fl", MESSVAR_NAME) ~ "above_canopy_front_left",
    grepl("oc_fr", MESSVAR_NAME) ~ "above_canopy_front_right",
    grepl("oc_bl", MESSVAR_NAME) ~ "above_canopy_back_left",
    grepl("oc_br", MESSVAR_NAME) ~ "above_canopy_back_right",
    grepl("uc_fl", MESSVAR_NAME) ~ "below_canopy_front_left",
    grepl("uc_fr", MESSVAR_NAME) ~ "below_canopy_front_right",
    grepl("uc_bl", MESSVAR_NAME) ~ "below_canopy_back_left",
    grepl("uc_br", MESSVAR_NAME) ~ "below_canopy_back_right",
    grepl("ic_fl", MESSVAR_NAME) ~ "in_canopy_front_left",
    grepl("ic_fr", MESSVAR_NAME) ~ "in_canopy_front_right",
    grepl("ic_bl", MESSVAR_NAME) ~ "in_canopy_back_left",
    grepl("ic_br", MESSVAR_NAME) ~ "in_canopy_back_right",
    grepl("ce", MESSVAR_NAME) ~ "center",
    TRUE ~ MESSVAR_NAME
  ))
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("above_canopy_front_left", "above_canopy_front_right",
                                                         "above_canopy_back_left", "above_canopy_back_right",
                                                         "in_canopy_front_left", "in_canopy_front_right",
                                                         "in_canopy_back_left", "in_canopy_back_right",
                                                         "center",
                                                         "below_canopy_front_left", "below_canopy_front_right",
                                                         "below_canopy_back_left", "below_canopy_back_right"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_ABS, 0.25)
Q3 <- quantile(filtered_data$DELTA_ABS, 0.75)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/AboveBelowInd_abs_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 6, height = 4)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_ABS, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (kPa)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) #+  # Dynamically set y-axis limits to exclude outliers
#scale_fill_manual(values = c("above_canopy" = "deepskyblue", "in_canopy" = "lightgreen", "center" = "yellow", "below_canopy" = "brown"))  # Custom color palette

dev.off()




###### figures VPD spread


### figure front to back relative change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("ic_ce", MESSVAR_NAME) & !grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("control", TREATMENT) ~ "control",
    grepl("irrigation", TREATMENT) ~ "irrigation",
    grepl("roof", TREATMENT) ~ "roof",
    TRUE ~ MESSVAR_NAME
  ))
# set filtered_data where DELTA = 0 to NA (remove where values = REF)
filtered_data <- filtered_data %>%
  mutate(
    DELTA_ABS = ifelse(DELTA_ABS == 0, NA, DELTA_ABS),
    DELTA_REL = ifelse(DELTA_REL == 0, NA, DELTA_REL)
  )
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("control", "irrigation", "roof"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_REL, 0.25, na.rm=T)
Q3 <- quantile(filtered_data$DELTA_REL, 0.75, na.rm=T)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/VPDSpread_rel_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_REL, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (%)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("control" = "yellow", "irrigation" = "deepskyblue", "roof" = "red"))  # Custom color palette

dev.off()



### figure front to back absolute change

# Filter the data for rows where "ic" appears in MESSVAR_NAME and "vpd" appears in TREATMENT
filtered_data <- meteo.df %>%
  filter(grepl("ic_ce", MESSVAR_NAME) & !grepl("vpd", TREATMENT))
# Create a new variable that combines 'bl' with 'br' and 'fl' with 'fr' in MESSVAR_NAME
filtered_data <- filtered_data %>%
  mutate(MESSVAR_NAME_COMBINED = case_when(
    grepl("control", TREATMENT) ~ "control",
    grepl("irrigation", TREATMENT) ~ "irrigation",
    grepl("roof", TREATMENT) ~ "roof",
    TRUE ~ MESSVAR_NAME
  ))
# set filtered_data where DELTA = 0 to NA (remove where values = REF)
filtered_data <- filtered_data %>%
  mutate(
    DELTA_ABS = ifelse(DELTA_ABS == 0, NA, DELTA_ABS),
    DELTA_REL = ifelse(DELTA_REL == 0, NA, DELTA_REL)
  )
# Define the order of the factor levels (front, center, back)
filtered_data$MESSVAR_NAME_COMBINED <- factor(filtered_data$MESSVAR_NAME_COMBINED, 
                                              levels = c("control", "irrigation", "roof"))
# Compute the limits for the y-axis to exclude outliers (using IQR)
Q1 <- quantile(filtered_data$DELTA_ABS, 0.25, na.rm=T)
Q3 <- quantile(filtered_data$DELTA_ABS, 0.75, na.rm=T)
IQR_value <- Q3 - Q1
lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/VPDSpread_abs_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 3, height = 3)

# Create the boxplot using ggplot2 with dynamically adjusted y-axis limits
ggplot(filtered_data, aes(x = MESSVAR_NAME_COMBINED, y = DELTA_ABS, fill = MESSVAR_NAME_COMBINED)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "black") +  # Add means as red points
  labs(y = expression(Delta~ " VPD (kPa)")) +  # Change y-axis label and remove title
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.title.x = element_blank(),  # Remove x-axis label
    legend.position = "none",  # Remove legend
    plot.title = element_blank()  # Remove plot title
  ) +
  scale_y_continuous(limits = c(lower_limit, upper_limit)) +  # Dynamically set y-axis limits to exclude outliers
  scale_fill_manual(values = c("control" = "yellow", "irrigation" = "deepskyblue", "roof" = "red"))  # Custom color palette

dev.off()




### figure operation time and mean VPD


# get hourly mean reference VPD

unique_median_ref <- meteo.df_complete %>%
  distinct(MESSTIME, MEDIAN_REF)
# Step 1: Extract the hour from MESSTIME (ignore date)
unique_median_ref <- unique_median_ref %>%
  mutate(hour = as.numeric(format(MESSTIME, "%H")))  # Extract hour (0-23)
# Step 2: Calculate the hourly mean of MEDIAN_REF
hourly_mean_ref <- unique_median_ref %>%
  group_by(hour) %>%
  summarise(mean_MEDIAN_REF = mean(MEDIAN_REF, na.rm = TRUE))
# Step 3: Create a data frame for all hours from 0 to 23
all_hours <- data.frame(hour = 0:23)
# Step 4: Merge the hourly mean data with the all_hours data frame (to ensure all hours 00:00 to 23:00 are included)
hourly_mean_ref_complete <- all_hours %>%
  left_join(hourly_mean_ref, by = "hour") %>%
  mutate(mean_MEDIAN_REF = ifelse(is.na(mean_MEDIAN_REF), 0, mean_MEDIAN_REF))  # Fill missing hours with 0


# get hourly mean manipulated VPD

unique_median_ref <- meteo.df_complete %>%
  distinct(MESSTIME, MEDIAN_VPD)
# Step 1: Extract the hour from MESSTIME (ignore date)
unique_median_ref <- unique_median_ref %>%
  mutate(hour = as.numeric(format(MESSTIME, "%H")))  # Extract hour (0-23)
# Step 2: Calculate the hourly mean of MEDIAN_REF
hourly_mean_ref <- unique_median_ref %>%
  group_by(hour) %>%
  summarise(mean_MEDIAN_REF = mean(MEDIAN_VPD, na.rm = TRUE))
# Step 3: Create a data frame for all hours from 0 to 23
all_hours <- data.frame(hour = 0:23)
# Step 4: Merge the hourly mean data with the all_hours data frame (to ensure all hours 00:00 to 23:00 are included)
hourly_mean_vpd_complete <- all_hours %>%
  left_join(hourly_mean_ref, by = "hour") %>%
  mutate(mean_MEDIAN_REF = ifelse(is.na(mean_MEDIAN_REF), 0, mean_MEDIAN_REF))  # Fill missing hours with 0



# plot

# Step 1: Get unique MESSTIME values (assuming meteo.df exists)
unique_meteo.df <- meteo.df %>%
  distinct(MESSTIME)
# Step 2: Extract hour from MESSTIME (ignore minute and second)
unique_meteo.df <- unique_meteo.df %>%
  mutate(hour = format(MESSTIME, "%H"))  # Extract hour
# Step 3: Count the frequency of each hour
hour_counts <- unique_meteo.df %>%
  count(hour) %>%
  mutate(hour = as.numeric(hour))  # Convert hour to numeric for correct sorting
# Step 4: Create a complete set of all hours (from 00 to 23)
all_hours <- data.frame(hour = 0:23)
# Step 5: Merge the counts with all possible hours (fill in missing hours with 0)
hour_counts_complete <- all_hours %>%
  left_join(hour_counts, by = "hour") %>%
  mutate(n = ifelse(is.na(n), 0, n))  # Replace NA values with 0
# Step 6: Calculate percentage of total hours
total <- length(seq(from = start_date, to = end_date, by = "days"))*60
hour_counts_complete <- hour_counts_complete %>%
  mutate(percentage = (n / total) * 100)

# save figure
pdf(paste0(getwd(), "/figures_", Sys.Date(), "/OpTimeMeanVPD_", format(start_date, "%Y-%m%-%d"), "_", format(end_date, "%Y-%m%-%d"), ".pdf"), width = 4.7, height = 4)

# Plot with secondary y-axis scaling and colored axis labels
ggplot(hour_counts_complete, aes(x = hour, y = percentage)) +
  geom_line(color = "blue", size = 1) +  # Frequency line plot
  geom_line(data = hourly_mean_ref_complete, 
            aes(x = hour, y = mean_MEDIAN_REF * 30),  # Scale mean_MEDIAN_REF by * 30
            color = "red", size = 1, linetype = "solid") +  # Mean MEDIAN_REF line plot
  geom_line(data = hourly_mean_vpd_complete, 
            aes(x = hour, y = mean_MEDIAN_REF * 30),  # Scale mean_MEDIAN_REF by * 30
            color = "red", size = 1, linetype = "dashed") +  # Mean MEDIAN_REF line plot
  labs(
    y = "Operation time VPD-manipulation (%)",  # Primary y-axis label
    x = "Hour (UTC)"  # Remove the x-axis label
  ) +
  scale_x_continuous(breaks = 0:23, labels = sprintf("%02d:00", 0:23), expand = c(0, 0)) +  # Hourly ticks on x-axis
  scale_y_continuous(
    name = "Operation time VPD-manipulation (%)",
    sec.axis = sec_axis(~ . / 30, name = "Mean VPD (kPa)")  # Scale secondary y-axis and update its label
  ) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(color = "blue"),  # Make primary y-axis labels blue
    axis.text.y.right = element_text(color = "red"),  # Make secondary y-axis labels red
    axis.title.y = element_text(color = "blue", size = 12),  # Primary y-axis title in blue
    axis.title.y.right = element_text(color = "red", size = 12),  # Secondary y-axis title in red
    panel.grid.major.x = element_line(color = "lightgray", size = 0.3, linetype = "solid"),  # Vertical gridlines at each hour
    panel.grid.minor.x = element_blank(),  # Remove minor vertical gridlines
    panel.grid.major.y = element_line(color = "lightgray", size = 0.3, linetype = "solid")  # Optional: Add horizontal gridlines
  )

dev.off()




### print important numbers

print("###")
print(paste0("Overall operation time from ", start_date, " to ", end_date, ": ", round(sum(meteo.df_complete$OPERATION_STATUS)/nrow(meteo.df_complete)*100,1), " %"))
print("###")
print(paste0("Delta-VPD from ", start_date, " to ", end_date, ": ", round(mean(meteo.df_complete$MEDIAN_VPD)- mean(meteo.df_complete$MEDIAN_REF), 2), " kPa"))
print("###")
print(paste0("Delta-VPD from ", start_date, " to ", end_date, ": ", round((mean(meteo.df_complete$MEDIAN_VPD)- mean(meteo.df_complete$MEDIAN_REF))/mean(meteo.df_complete$MEDIAN_VPD)*100, 1), " %"))
print("###")
print(paste0("Delta-VPD during operation time from ", start_date, " to ", end_date, ": ", round(meteo.df_complete %>%
                                                                                                  filter(OPERATION_STATUS == 1) %>%
                                                                                                  summarise(mean_vpd = mean(MEDIAN_VPD)) %>%
                                                                                                  pull(mean_vpd)
                                                                                                - meteo.df_complete %>%
                                                                                                  filter(OPERATION_STATUS == 1) %>%
                                                                                                  summarise(mean_vpd = mean(MEDIAN_REF)) %>%
                                                                                                  pull(mean_vpd)
                                                                                                , 2), " kPa"))
print("###")
print(paste0("Delta-VPD during operation time from ", start_date, " to ", end_date, ": ", round((meteo.df_complete %>%
                                                                                                   filter(OPERATION_STATUS == 1) %>%
                                                                                                   summarise(mean_vpd = mean(MEDIAN_VPD)) %>%
                                                                                                   pull(mean_vpd)
                                                                                                 - meteo.df_complete %>%
                                                                                                   filter(OPERATION_STATUS == 1) %>%
                                                                                                   summarise(mean_vpd = mean(MEDIAN_REF)) %>%
                                                                                                   pull(mean_vpd))
                                                                                                / meteo.df_complete %>%
                                                                                                  filter(OPERATION_STATUS == 1) %>%
                                                                                                  summarise(mean_vpd = mean(MEDIAN_REF)) %>%
                                                                                                  pull(mean_vpd) *100
                                                                                                , 1), " %"))

