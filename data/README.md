Chicago.RData: Chicago crime dataset consists of crime counts reported in the city of Chicago, ranging from January
1st, 2001 to December 11th, 2017. The follwings are variables in the RData.

* ltns: An order-3 tensor with entries representing the log counts of crimes from 24 hours, 77 community
areas, and 32 crime types.
* hour_map: labels for 24 hours.
* area_map: labels for 77 community areas. 
* crimetype_map: labels for 32 crime types.

Chicago_result.R: Code to generate the results in the main paper, which include Figure 6-7, Table 4, and Table S3.
Chicago_data_processing: Code for data cleaning from the source http://frostt.io/tensors/chicago-crime/

