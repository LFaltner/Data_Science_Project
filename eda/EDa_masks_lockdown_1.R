# clear brain

rm(list= ls())

# import of libraries
library("ggplot2")
library("dplyr")
library("ggfortify")
library('tidyverse')
library('GGally')

# load data

setwd("~/BSc_UZH/UZH_22FS/ESC403_Project/wrangled_data")

cases <- read.csv('cases_CH.csv', stringsAsFactors = T)
lockdown <- read.csv('lockdown_CH.csv', stringsAsFactors = T)
masks <- read.csv('masks_CH.csv', stringsAsFactors = T)



# rough idea

glimpse(cases)
glimpse(masks)
glimpse(lockdown)

cases <- slice(cases, 1:766)
# distribution of cases

ggplot(data = cases) +
  geom_histogram(mapping = aes(x = new_cases), bins = 100)

# cases against time line

ggplot(data = cases, mapping = aes(x = seq(as.Date('2020-02-25'), as.Date('2022-03-31'), by = '1 day'), y = new_cases_smoothed, color = masks$facial_coverings)) +
  geom_line()

ggplot(data = cases, mapping = aes(x = seq(as.Date('2020-02-25'), as.Date('2022-03-31'), by = '1 day'), y = new_cases_smoothed, color = lockdown$stay_home_requirements)) +
  geom_line()

  

