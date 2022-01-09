library(vegan)
library(LaplacesDemon)

raw = read_xlsx("H:/JenaSP6_2021/id034_nematodes_240-samples_JenaExperiment_DeltaBEF_July2021_Marcel-Ciobanu__240 samples Jena exp SP6 Nematodes Amyntas2021FINAL.xlsx",
                sheet = "Synthesis",
                range = "A5:BM245")
raw[is.na(raw)] = 0

# arrange by sample
data = raw %>% arrange(Sample) 

# replace underscore with D to match the Jexis naming scheme 
data$Sample = gsub("_", "D", data$Sample) 

# break that into block, plot, treatment
this = str_split(data$Sample, "A|D", simplify = T)

data.2 = data %>% add_column(block = this[,1],
                             plot = this[,2],
                             treatment = this[,3],
                             .after = "Sample")

# the unlabeled sample is B1A12D3
table(data.2$block,data.2$plot)
View(filter(data.2, block == "B1" & plot == 12))
data.2[240,1:4] = as.list(c("B1A12D3","B1","12","3"))

# arrange again
data.2 = data.2 %>% arrange(Sample)

# get nematode abundances & soil info
gs4_deauth() #does this work?
abun = read_sheet("https://docs.google.com/spreadsheets/d/1YMQmyhLYfr86CcmwpLwRkLQPxwy4Yk2oyrPKXO8Cf0w/edit#gid=0",
                  sheet = "abundance")

data.3 = data.2[,5:68]
rownames(data.3) = data.2$Sample

View(data.3[,12])

data.4 = as.matrix(vegan::decostand(data.3,"total",1))

extra = vector(mode="list",100)
for (i in 1:100) {
  extra[[i]] = matrix(0, nrow(data.4),ncol(data.4))
  temp.1 = vector(mode = "list", nrow(data.4))
  for (j in 1:nrow(data.4)) { 
  temp.1[[j]] = as.data.frame(table(sample(colnames(data.4),
                                         as.numeric(abun[j,4]),
                                         replace=TRUE, 
                                         prob=data.4[j,]))) 
  names(temp.1[[j]]) = c("taxon", rownames(data.4)[j])
  }
  temp.2 = temp.1 %>% reduce(full_join, by = "taxon")
  rownames(temp.2) = temp.2[,1] ; temp.2 = t(temp.2[,-1]) ; temp.2[is.na(temp.2)] = 0
  extra[[i]] = temp.2
}

extrapolate.robust = function(df, counts, n = 100) {
  # create a list of length equal to the number of "versions"
  extra = vector(mode="list", n)
  for (i in 1:n) {
    # each element of the list is an empty matrix
    extra[[i]] = matrix(0, nrow(df),ncol(df))
    # another list to hold info on the random draws for each site
    temp.1 = vector(mode = "list", nrow(df))
    for (j in 1:nrow(df)) { 
      # we draw randomly, with replacement, from the available taxa
      # for each site the number of draws is equal to the number of counted nematodes
      # the probability of a taxon being drawn is the frequency of that taxon 
      # in the identified nematodes for that site
      temp.1[[j]] = as.data.frame(table(sample(colnames(df), 
                                               as.numeric(counts), 
                                               replace=TRUE, 
                                               prob=df[j,]))) 
      # internally, the above produces one dataframe per site 
      # where the first column is (only) the drawn taxa and the second one 
      # the numbers of those taxa
      # we name them "taxon" and sample name respectively
      names(temp.1[[j]]) = c("taxon", rownames(df)[j])
    }
    # now we join the dataframes of all sites into one
    temp.2 = temp.1 %>% reduce(full_join, by = "taxon")
    # NAs are artifacts of joining, placed where a taxon is missing from a site
    temp.2[is.na(temp.2)] = 0
    
    rownames(temp.2) = temp.2[,1] # we name rows after taxa,
    temp.2 = t(temp.2[,-1])  # drop the column that contained the taxa names and transpose
    # the resulting sites by species matrix goes to the list of "versions"
    extra[[i]] = temp.2
  }
  
  return(extra)
}

well = extrapolate.robust(data.4, abun[j,4], 100)



library(tidyverse)
content <- read_html("http://nemaplex.ucdavis.edu/Ecology/bioindicators.htm")
install.packages("xml2")
library(xml2)
content <- read_html("http://nemaplex.ucdavis.edu/Ecology/bioindicators.htm")
View(content)
View(content$node)
View(content$html)
View(content$html$body$table)
View(content[[1]][2])
View(content[[1]][[2]])
tables <- content %>% html_table(fill = TRUE)
library(tidyverse)
tables <- content %>% html_table(fill = TRUE)
install.packages("rvest")
library(rvest)
content <- read_html("http://nemaplex.ucdavis.edu/Ecology/bioindicators.htm")
tables <- content %>% html_table(fill = TRUE)
View(tables)
View(tables[[1]])
View(tables[[2]])
first = tables[[2]][-1,1:4]
View(first)
first = tables[[2]][-1,1:3]
View(first)
View(first)
?rename
library(tidyverse)
?rename
first = tables[[2]][-1,1:3] %>% rename(X1="Family")
first = tables[[2]][-1,1:3] %>% rename(Family=X1)
View(first)
first = tables[[2]][-1,1:3] %>% rename(c(Family,`Feeding Habit`)=c(X1,X2))
colnames(first) = first[1,]
View(first)
first = first[-1,]
View(first)
first = tables[[2]][-1,1:3]
colnames(first) = first[1,]
View(first)
first = first[-1,]
second = tables[[3]][-1,1:3]
colnames(second) = second[1,]
second = second[-1,]
View(second)
third = tables[[4]][-1,1:3]
colnames(third) = third[1,]
third = third[-1,]
View(third)
?rbind
all = rbind(first, second, third)




string = "AchromadoridaeActinolaimidaeAlaimidaeAlirhabditidaeAlloionematidaeAmphidelidaeAnatonchidaeAnguinidaeAphanolaimidaeAphelenchidaeAphelenchoididaeAporcelaimidaeAulolaimidaeAulolaimoididaeBastianiidaeBathyodontidaeBelondiridaeBicirronematidaeBunonematidaeCamacolaimidaeCampydoridaeCephalobidaeChambersiellidaeChoanolaimidaeChromadoridaeChronogastridaeCobbonchidaeCrateronematidaeCriconematidaeCyatholaimidaeDesmodoridaeDesmoscolecidaeDiphtherophoridaeDiplogasteroididaeDiplogastridaeDiplopeltidaeDiploscapteridaeDolichodoridaeDorylaimellidaeDorylaimidaeElaphonematidaeEncholaimidaeEthmolaimidaeHemicycliophoridaeHeteroderidaeHoplolaimidaeHypodontolaimidaeIotonchidaeIotonchiidaeIronidaeIsolaimiidaeLeptolaimidaeLeptonchidaeLinhomoeidaeLongidoridaeMerliniidaeMicrolaimidaeMonhysteridaeMononchidaeMononchulidaeMydonomidaeMylonchulidaeMyolaimidaeNeodiplogastridaeNeotylenchidaeNordiidaeNygellidaeNygolaimidaeOdontolaimidaeOdontopharyngidaeOhrididaeOnchulidaeOsstellidaePanagrolaimidaeParatylenchidaeParaxonchiidaePlectidaePratylenchidaePrismatolaimidaePterygorhabditidaeQudsianematidaeRhabditidaeRhabdolaimidaeRhadinematidaeSphaerolaimidaeSphaeronematidaeSphaerulariidaeSteinernematidaeSwangeriidaeTelotylenchidaeTeratocephalidaeThoracostomopsidaeThornenematidaeThorniidaeTobrilidaeTrichodoridaeTripylidaeTrischistomatidaeTyencholaimidaeTylenchidaeTylencholaimellidaeTylencholaimidaeTylenchulidaeTylopharyngidaeXyalidae"
taxa = strsplit(string, "(?<=[a-z])(?=[A-Z])", perl = TRUE)[[1]]






library(rvest)
library(rJava)
library(RSelenium)

driver <- rsDriver(browser = "firefox")
remDr<-driver[["client"]]

remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/FamilyParmsQuery.aspx") 

opt1 <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList2"]')
opt1$clickElement()

opt2 <- remDr$findElement(using = 'xpath', "//*/option[@value = 'Achromadoridae']")
opt2$clickElement()

opt3 <- remDr$findElement(using = 'xpath', value= '//*[@id="Button1"]')
opt3$clickElement()


#get tables
find.tables <- remDr$findElement(using = 'xpath', value= '//*[@id="DetailsView5"]')

tables = tb$getPageSource()[[1]] %>% read_html() %>% html_table(fill = TRUE)
dada = dada[[4]]


View(rbind(tables[[4]], tables[[6]], tables[[6]]))


dadada = rbind(tables[[4]],tables[[5]],tables[[6]]) %>% 
         setNames(dadada[1,]) %>% dplyr::slice(-1)

 e      
element = paste(c(substr(element, 1, 21), 
                  substr(element, 36,38)), 
                collapse = strings[2])


library(tidyverse)
library(rvest)
library(rJava)
library(RSelenium)


string = "AchromadoridaeActinolaimidaeAlaimidaeAlirhabditidaeAlloionematidaeAmphidelidaeAnatonchidaeAnguinidaeAphanolaimidaeAphelenchidaeAphelenchoididaeAporcelaimidaeAulolaimidaeAulolaimoididaeBastianiidaeBathyodontidaeBelondiridaeBicirronematidaeBunonematidaeCamacolaimidaeCampydoridaeCephalobidaeChambersiellidaeChoanolaimidaeChromadoridaeChronogastridaeCobbonchidaeCrateronematidaeCriconematidaeCyatholaimidaeDesmodoridaeDesmoscolecidaeDiphtherophoridaeDiplogasteroididaeDiplogastridaeDiplopeltidaeDiploscapteridaeDolichodoridaeDorylaimellidaeDorylaimidaeElaphonematidaeEncholaimidaeEthmolaimidaeHemicycliophoridaeHeteroderidaeHoplolaimidaeHypodontolaimidaeIotonchidaeIotonchiidaeIronidaeIsolaimiidaeLeptolaimidaeLeptonchidaeLinhomoeidaeLongidoridaeMerliniidaeMicrolaimidaeMonhysteridaeMononchidaeMononchulidaeMydonomidaeMylonchulidaeMyolaimidaeNeodiplogastridaeNeotylenchidaeNordiidaeNygellidaeNygolaimidaeOdontolaimidaeOdontopharyngidaeOhrididaeOnchulidaeOsstellidaePanagrolaimidaeParatylenchidaeParaxonchiidaePlectidaePratylenchidaePrismatolaimidaePterygorhabditidaeQudsianematidaeRhabditidaeRhabdolaimidaeRhadinematidaeSphaerolaimidaeSphaeronematidaeSphaerulariidaeSteinernematidaeSwangeriidaeTelotylenchidaeTeratocephalidaeThoracostomopsidaeThornenematidaeThorniidaeTobrilidaeTrichodoridaeTripylidaeTrischistomatidaeTyencholaimidaeTylenchidaeTylencholaimellidaeTylencholaimidaeTylenchulidaeTylopharyngidaeXyalidae"
taxa = strsplit(string, "(?<=[a-z])(?=[A-Z])", perl = TRUE)[[1]]

element = "//*/option[@value = 'Achromadoridae']"

driver <- rsDriver(browser = "firefox")
remDr<-driver[["client"]]

dadada = vector(mode = "list", length(taxa))
#names(dadada) = taxa

for (taxon in 1:length(taxa)) {
  remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/FamilyParmsQuery.aspx") 
  
  opt1 <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList2"]')
  opt1$clickElement()
  
  opt2 <- remDr$findElement(using = 'xpath', paste(c(substr(element, 1, 21), 
                                                     substr(element, 36,38)), 
                                                   collapse = taxa[taxon]))
  opt2$clickElement()
  
  opt3 <- remDr$findElement(using = 'xpath', value= '//*[@id="Button1"]')
  opt3$clickElement()
  
  
  #get tables
  find.tables <- remDr$findElement(using = 'xpath', value= '//*[@id="DetailsView5"]')
  
  tables = find.tables$getPageSource()[[1]] %>% read_html() %>% html_table(fill = TRUE)
  #tables = tables[[4]]
  
  dadada[[taxon]] = rbind(tables[[4]],tables[[5]],tables[[6]]) #%>% 
                    #setNames(dadada[[taxon]][1,]) #%>% dplyr::slice(-1)
  
  dadada[[taxon]] = dadada[[taxon]] %>% 
                    setNames(dadada[[taxon]][1,]) %>% 
                    dplyr::slice(-1)
  
  
}

da.da.da = dadada %>% reduce(full_join, by = "Family")
da.da.da = as.data.frame(da.da.da)
rownames(da.da.da) = da.da.da[,1]
da.da.da = da.da.da[,-1]
ad.ad.ad = as.data.frame(t(da.da.da))

# Save a single object to a file
saveRDS(ad.ad.ad, "nemaplex_family.rds")
# Restore it under a different name
my_data <- readRDS("nemaplex_family.rds")




remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/GenusParmsQuery.aspx")
webElem <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList3"]')
genera <- webElem$getElementText()
taxa = strsplit(genera[[1]], "\n", perl = TRUE)[[1]]

dadada = vector(mode = "list", length(taxa))
#names(dadada) = taxa

for (taxon in 1:length(taxa)) {
  remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/GenusParmsQuery.aspx") 
  
  opt1 <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList3"]')
  opt1$clickElement()
  
  opt2 <- remDr$findElement(using = 'xpath', paste(c(substr(element, 1, 21), 
                                                     substr(element, 36,38)), 
                                                   collapse = taxa[taxon]))
  opt2$clickElement()
  
  opt3 <- remDr$findElement(using = 'xpath', value= '//*[@id="Button1"]')
  opt3$clickElement()
  
  
  #get tables
  find.tables <- remDr$findElement(using = 'xpath', value= '//*[@id="DetailsView2"]')
  
  tables = find.tables$getPageSource()[[1]] %>% read_html() %>% html_table(fill = TRUE)
  #tables = tables[[4]]
  
  dadada[[taxon]] = rbind(tables[[3]],tables[[4]],tables[[5]]) #%>% 
  #setNames(dadada[[taxon]][1,]) #%>% dplyr::slice(-1)
  
  dadada[[taxon]] = dadada[[taxon]] %>% 
    setNames(dadada[[taxon]][1,]) %>% 
    dplyr::slice(-1)
  
  
}

da.da.da = dadada %>% reduce(full_join, by = "Genus")
da.da.da = as.data.frame(da.da.da)
rownames(da.da.da) = da.da.da[,1]
da.da.da = da.da.da[,-1]
ad.ad.ad = as.data.frame(t(da.da.da))

# Save a single object to a file
saveRDS(ad.ad.ad, "nemaplex_genus.rds")




getttt = query_nemaplex(c("Achromadoridae","Actinolaimidae"))



options(width = 80)
n <- 100
for (ii in 1:n) {
  extra <- nchar('||100%')
  width <- options()$width
  step <- round(ii / n * (width - extra))
  text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
                  strrep(' ', width - step - extra), round(ii / n * 100))
  cat(text)
  Sys.sleep(0.05)
  cat(if (ii == n) '\n' else '\014')
}


getttt$StDevMass = getttt$StderrMass * sqrt(getttt$N)
108.3382
27.4822
2.829466e+02
hist(rlnormTrunc(10000, 
                 meanlog = 108.3382, 
                 sdlog = 2.829466e+02, 
                 min = qlnorm(.05, 
                              meanlog = 108.3382, 
                              sdlog = 2.829466e+02), 
                 max = qlnorm(.95, 
                              meanlog = 108.3382, 
                              sdlog = 2.829466e+02)),
     breaks = 1000)

hist(rlnorm(1000, 
            meanlog = 108.3382, 
            sdlog = 2.829466e+02),
     breaks = 100)


bla = rlnorm(10000, 
       meanlog = 0.3913, 
       sdlog = 4.242641e-03
)
hist(log(bla),breaks = 1000)

hist( rtruncnorm( 50000, 
                  a=qnorm(.05, log(0.5493   ), 
                               log(1.668325 )), 
                  b=qnorm(.95, log(0.5493   ), 
                               log(1.668325 )), 
                        mean = log(0.5493   ), 
                          sd = log(1.668325 )), 
      breaks = 1000)


m <- 7
s <- 75
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
print(paste("location:", location))
print(paste("shape:", shape))
draws3 <- rlnorm(n=1000000, location, shape)

hist(draws3, 
     breaks = 1000)



m <- 26.4783
s <- 3.821168e+01
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))
hist(rlnormTrunc(10000, 
                 meanlog = location, 
                 sdlog = shape, 
                 min = qlnorm(.05, 
                              meanlog = location, 
                              sdlog = shape), 
                 max = qlnorm(.95, 
                              meanlog = location, 
                              sdlog = shape)),
     breaks = 1000)
