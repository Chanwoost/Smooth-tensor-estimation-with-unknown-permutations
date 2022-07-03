data=read.table("rawdata/chicago-crime-comm.tns")
date_map=read.table("mode-1-date.map")
hour_map=read.csv("mode-2-hour.map",head = F)
area_map=read.csv("mode-3-communityarea.map",head = F)
crimetype_map=read.csv("mode-4-crimetype.map",head = F)
# reshape dataframe to tensor
names(data) = c("date","hour","area","type","count")
ndata = aggregate(data[5],data[2:4],sum)


library(reshape2)
tns = acast(ndata,hour~area~type,value.var = "count")
tns[is.na(tns)] = 0
ltns = log(1+tns)

save(tns,ltns,hour_map,area_map,crimetype_map,file = "Chicago.RData")