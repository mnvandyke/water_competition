rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(moments)

precip <- read_csv("./data/2019_rainfall.csv")
precip$month <- factor(precip$month, levels = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"))


ggplot(precip) +
  geom_point(aes(x = month, y = cumulative_dry_plots, col = "Dry Plots", size = 0.25)) +
  geom_point(aes(x = month, y = cumulative_wet_plots, col = "Wet Plots", size = 0.25)) +
  scale_color_manual(values=c("Wet Plots" ="#00BFC4", "Dry Plots" = "#F8766D")) +
  ylab("Cumulative precipitation (cm)") +
  ggtitle("Precipitation Treatments")
ggsave("./figures/precip_treatments.pdf")

fig_mtn_rain <- read_csv("./data/Fig_mtn_rain_data.csv")  
ggplot(fig_mtn_rain,aes(x = total_cm )) +
  geom_histogram(binwidth = 2,  fill="deepskyblue3") +
  geom_vline(xintercept = mean(fig_mtn_rain$total_cm), size = 1) +
  #geom_vline(xintercept = median(fig_mtn_rain$total_cm), color = "green") +
  geom_vline(xintercept = 52.18, color= "red", size = 1) +
  geom_vline(xintercept = 64.94, color = "blue", size = 1) +
  ylab("Frequency of years") +
  xlab("Rain Fall (cm)") +
  ggtitle("Yearly Rainfall totals 1960-2018")
ggsave("./figures/rainfall_hist.pdf")

df2<- ("sn"= (c(2,3.5)))
mean(fig_mtn_rain$total_cm)
ggplot() +
  geom_point(x = 2, y = 5) +
  geom_abline( intercept = 0, slope = 1) +
  ylim(0, 10) +
  xlim(0, 1)

ggplot(fig_mtn_rain,aes(x = year_start, y = total_cm ))+
  geom_point() +
  geom_line( size = 1) +
  ylab("Total rainfall (cm)") +
  xlab("Year") +
  ylim(0, 135) +
  theme(axis.title  = element_text(vjust=0.5, size=16))
ggsave("./figures/fig_rainfall.pdf", height = 4, width = 7)

# Which months have the most variation in rainfall?
precip <- read.csv("./data/daily_rainfall_data_figmtn.csv")

months <- precip %>%
  group_by(year) %>%
  group_by(month) %>%
  summarise(sum(daily.rain))

jan_precip <- precip %>%
  filter(month == 1)
jan_years<- jan_precip %>%
  group_by(year) %>%
  summarise(sum(daily.rain))

mean(jan_years$`sum(daily.rain)`)
sd(jan_years$`sum(daily.rain)`)
kurtosis(jan_years$`sum(daily.rain)`)

feb_precip <- precip %>%
  filter(month == 2)
feb_years<- feb_precip %>%
  group_by(year) %>%
  summarise(sum(daily.rain))

mean(feb_years$`sum(daily.rain)`)
sd(feb_years$`sum(daily.rain)`)
kurtosis(feb_years$`sum(daily.rain)`)

mar_precip <- precip %>%
  filter(month == 3)
mar_years<- mar_precip %>%
  group_by(year) %>%
  summarise(sum(daily.rain))

mean(mar_years$`sum(daily.rain)`)
sd(mar_years$`sum(daily.rain)`)
kurtosis(mar_years$`sum(daily.rain)`)

apr_precip <- precip %>%
  filter(month == 4)
apr_years<- apr_precip %>%
  group_by(year) %>%
  summarise(sum(daily.rain))

mean(apr_years$`sum(daily.rain)`)

sd(apr_years$`sum(daily.rain)`)
kurtosis(apr_years$`sum(daily.rain)`)

may_precip <- precip %>%
  filter(month == 5)
may_years<- may_precip %>%
  group_by(year) %>%
  summarise(sum(daily.rain))

mean(may_years$`sum(daily.rain)`)

sd(may_years$`sum(daily.rain)`)
kurtosis(may_years$`sum(daily.rain)`)


#Categorizing rainfall years into above and below average

fig_mtn_rain$rain_type <- ifelse(fig_mtn_rain$total_cm>= mean(fig_mtn_rain$total_cm), "above", "below")

