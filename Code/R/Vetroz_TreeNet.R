# TreeNet data in Vetroz

# read in timeseries and metadata
TN = read_csv("TreeNet_dendro/tn_timeseries.csv")
meta = read_csv("TreeNet_dendro/tn_metadata.csv")

comp_plots = function(TN, id1, id2) {
  
  TN1 = filter(TN, series_id == id1)
  TN2 = filter(TN, series_id == id2)
  
  TN1_clip = TN1 %>% filter(ts >= TN2$ts[1]) %>% 
    rename(value1 = value) %>% select(ts, value1)
  TN2 = left_join(TN2, TN1_clip)
  
  p = ggplot(TN2, aes(ts, value, color=as.character(id1)))+geom_point(size=.5)+
    geom_point(aes(ts, value1, color=as.character(id2)), size=.5)
  
  ggsave(sprintf("Series%d_%dcomp.png", id1, id2), plot=p, device="jpeg")
  
}

# first, compare duplicates
# series 1357 (sensor 15417) is worst
TN1357 = filter(TN, series_id == 1357) # longer record with unreasonable values
TN1812 = filter(TN, series_id == 1812) # shorter record with reasonable values

ggplot(TN1357, aes(ts, value))+geom_point(size=.5)
ggplot(TN1812, aes(ts, value))+geom_point(size=.5)

TN1357_clip = TN1357 %>% filter(ts >= TN1812$ts[1]) %>% 
  rename(value1357 = value) %>% select(ts, value1357)
TN1812 = left_join(TN1812, TN1357_clip)

TN1812$diff = TN1812$value - TN1812$value1357
TN1812$value1357_fix = TN1812$value1357 + mean(TN1812$diff, na.rm=T)

ggplot(TN1812, aes(ts, value, color="1812"))+geom_point(size=.5)+
  geom_point(aes(ts, value1357_fix, color="1357"), size=.5)+
  geom_point(aes(ts, diff, color="diff"), size=.5)

TN1357$value_fix = TN1357$value + TN1812$diff[1]
ggplot(TN1357, aes(ts, value_fix))+geom_point(size=.5)


# sensor 15418
TN1362 = filter(TN, series_id == 1362) # longer record
TN1813 = filter(TN, series_id == 1813) # shorter record

ggplot(TN1362, aes(ts, value))+geom_point(size=.5)
ggplot(TN1813, aes(ts, value))+geom_point(size=.5)

TN1362_clip = TN1362 %>% filter(ts >= TN1813$ts[1]) %>% 
  rename(value1362 = value) %>% select(ts, value1362)
TN1813 = left_join(TN1813, TN1362_clip)

TN1813$diff = TN1813$value1362 - TN1813$value

ggplot(TN1813, aes(ts, value, color="1813"))+geom_point(size=.5)+
  geom_point(aes(ts, value1362, color="1362"), size=.5)


# sensor 15419
TN1363 = filter(TN, series_id == 1363) # longer record
TN1814 = filter(TN, series_id == 1814) # shorter record

ggplot(TN1363, aes(ts, value))+geom_point(size=.5)
ggplot(TN1814, aes(ts, value))+geom_point(size=.5)

TN1363_clip = TN1363 %>% filter(ts >= TN1814$ts[1]) %>% 
  rename(value1363 = value) %>% select(ts, value1363)
TN1814 = left_join(TN1814, TN1363_clip)

#TN1814$diff = TN1813$value1363 - TN1814$value

ggplot(TN1814, aes(ts, value, color="1814"))+geom_point(size=.5)+
  geom_point(aes(ts, value1363, color="1363"), size=.5)


comp_plots(TN, 1363, 1814)
# sensor 15419
TN1363 = filter(TN, series_id == 1363) # longer record
TN1814 = filter(TN, series_id == 1814) # shorter record

ggplot(TN1363, aes(ts, value))+geom_point(size=.5)
ggplot(TN1814, aes(ts, value))+geom_point(size=.5)

TN1363_clip = TN1363 %>% filter(ts >= TN1814$ts[1]) %>% 
  rename(value1363 = value) %>% select(ts, value1363)
TN1814 = left_join(TN1814, TN1363_clip)

#TN1814$diff = TN1813$value1363 - TN1814$value

ggplot(TN1814, aes(ts, value, color="1814"))+geom_point(size=.5)+
  geom_point(aes(ts, value1363, color="1363"), size=.5)


comp_plots(TN, 1358, 1815)
# sensor 15422
TN1358 = filter(TN, series_id == 1358) # longer record
TN1815 = filter(TN, series_id == 1815) # shorter record

ggplot(TN1358, aes(ts, value))+geom_point(size=.5)
ggplot(TN1815, aes(ts, value))+geom_point(size=.5)

TN1358_clip = TN1358 %>% filter(ts >= TN1815$ts[1]) %>% 
  rename(value1358 = value) %>% select(ts, value1358)
TN1815 = left_join(TN1815, TN1358_clip)

#TN1814$diff = TN1813$value1363 - TN1814$value

ggplot(TN1815, aes(ts, value, color="1815"))+geom_point(size=.5)+
  geom_point(aes(ts, value1358, color="1358"), size=.5)


comp_plots(TN, 1359, 1816)


TN1364 = filter(TN, series_id == 1364)
ggplot(TN1364, aes(ts, value))+geom_point()