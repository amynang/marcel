# this one needs checking and tests 
extrapolate = function(df, counts) {
  extra = df %>% mutate(across(everything(), # for all columns
                               # we multiply by total abundance and divide by 100
                               ~ .*as.numeric(counts)/100))
  return(extra)
}



extrapolate.robust = function(df, counts, n = 100) {
  
  require(magrittr)
  
  if (!is.numeric(df)) {
    stop("df must be numeric")
  }
  if (sum(is.na(df))>0) {
    stop("df contains NAs")
  }
  
  # if the df is raw counts we get frequencies
  df.f = vegan::decostand(df,"total",1)
  
  # create a list of length equal to the number of "versions"
  extra = vector(mode="list", n)
  for (i in 1:n) {
    # each element of the list is an empty matrix
    extra[[i]] = matrix(0, nrow(df.f),ncol(df.f))
    # another list to hold info on the random draws for each site
    temp.1 = vector(mode = "list", nrow(df.f))
    for (j in 1:nrow(df.f)) { 
      # we draw randomly, with replacement, from the available taxa
      # for each site the number of draws is equal to the number of counted nematodes
      # the probability of a taxon being drawn is the frequency of that taxon 
      # in the identified nematodes for that site
      temp.1[[j]] = as.data.frame(table(sample(colnames(df.f), 
                                               as.numeric(counts), 
                                               replace=TRUE, 
                                               prob=df.f[j,]))) 
      # internally, the above produces one dataframe per site 
      # where the first column is (only) the drawn taxa and the second one 
      # the numbers of those taxa
      # we name them "taxon" and sample name respectively
      names(temp.1[[j]]) = c("taxon", rownames(df.f)[j])
    }
    # now we join the dataframes of all sites into one
    temp.2 = temp.1 %>% purrr::reduce(dplyr::full_join, by = "taxon")
    # NAs are artifacts of joining, placed where a taxon is missing from a site
    temp.2[is.na(temp.2)] = 0
    
    rownames(temp.2) = temp.2[,1] # we name rows after taxa,
    temp.2 = t(temp.2[,-1])  # drop the column that contained the taxa names and transpose
    # the resulting sites by species matrix goes to the list of "versions"
    extra[[i]] = temp.2
    
    # show loop progress
    cat('\014')
    cat(paste0(round((i/n)*100), '% completed'))
    Sys.sleep(.05)
    if (i == n) cat(': Done')
    
  }
  
  return(extra)
}


query_nemaplex = function(taxa, complete = FALSE) {
  require(tidyverse)
  require(rvest)
  require(rJava)
  require(RSelenium)
  
  # start the thingy
  driver <- rsDriver(browser = "firefox", verbose = FALSE
                     #this line makes a headless firefox :) 
                     #i.e. running in background; comment it out for debugging
                     ,extraCapabilities = list("moz:firefoxOptions" = list(args = list('--headless')))
                     )
  remDr<-driver[["client"]]
  
  
  ############################## Family level Query ##############################
  
  # which url
  remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/FamilyParmsQuery.aspx")
  # which element inside that url (on Firefox, you can find that by 
  # rigth click > Inspect and then patiently going through the page elements)
  webElem <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList2"]')
  # get the list of available families on the dropdown menu
  families <- webElem$getElementText()
  f.present = strsplit(families[[1]], "\n", perl = TRUE)[[1]]
  
  # filter your taxa to get only those that are in the above list
  taxa.s = taxa[taxa %in% f.present]
  
  if (length(taxa.s) > 0) { 
    
    # taxa names will be interposed between the ''
    element = "//*/option[@value = '']"
    
    family.list = vector(mode = "list", length(taxa.s))
    #names(dadada) = taxa
    
    for (taxon in 1:length(taxa.s)) {
      # which url
      remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/FamilyParmsQuery.aspx") 
      # which element
      opt1 <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList2"]')
      # click
      opt1$clickElement()
      # we place a Family between the ''...
      opt2 <- remDr$findElement(using = 'xpath', paste(c(substr(element, 1, 21), 
                                                         substr(element, 22,23)), 
                                                       collapse = taxa.s[taxon]))
      # ...then click
      opt2$clickElement()
      # we select the Process Query button and click
      opt3 <- remDr$findElement(using = 'xpath', value= '//*[@id="Button1"]')
      opt3$clickElement()
      
      # ta-dan! we are inside the Family Parameters page
      
      # I do not know why this works but it gives a list with all tables
      find.tables <- remDr$findElement(using = 'xpath', value= '//*[@id="DetailsView5"]')
      tables = find.tables$getPageSource()[[1]] %>% read_html() %>% html_table(fill = TRUE)
      
      # we select the tables with taxon parameters and combine them
      family.list[[taxon]] = rbind(tables[[4]],tables[[5]],tables[[6]]) 
      
      # name columns after 1st row ("Family","(family name)")
      family.list[[taxon]] = family.list[[taxon]] %>% 
        setNames(family.list[[taxon]][1,]) %>% 
        # then drop that row
        dplyr::slice(-1)
      
      cat('\014')
      cat(paste0(round((taxon/length(taxa.s))*100), '%'))
      #Sys.sleep(.05)
      if (taxon == length(taxa.s)) cat('- Families query complete
                                       ')
      
    }
    
    taxa.fam = family.list %>% reduce(full_join, by = "Family")
    taxa.fam = as.data.frame(taxa.fam)
    rownames(taxa.fam) = taxa.fam[,1]
    taxa.fam = taxa.fam[,-1, drop=FALSE]
    all.fam = as.data.frame(t(taxa.fam))
    all.fam[,5:28] = data.frame(lapply(all.fam[,5:28], as.numeric))
    all.fam = dplyr::rename(all.fam, N=nFamily)
    names(all.fam) = gsub("Famavg", "Avg", names(all.fam))
  }else{
    # in case your input contains no families
    all.fam = vector(mode = "list", length(taxa.s))
    print("No families in input")
  }
  
  ############################## Genus level Query ###############################
  
  # which url
  remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/GenusParmsQuery.aspx")
  # which element inside that url (on Firefox, you can find that by 
  # rigth click > Inspect and then patiently going through the page elements)
  webElem <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList3"]')
  # get the list of available genera on the dropdown menu
  genera <- webElem$getElementText()
  g.present = strsplit(genera[[1]], "\n", perl = TRUE)[[1]]
  
  # filter your taxa to get only those in the above list
  taxa.s = taxa[taxa %in% g.present]
  
  if (length(taxa.s) > 0) { 
    
    # taxa names will be interposed between the ''
    element = "//*/option[@value = '']"
    
    genus.list = vector(mode = "list", length(taxa.s))
    
    
    for (taxon in 1:length(taxa.s)) {
      # which url
      remDr$navigate("http://nemaplex.ucdavis.edu/Ecology/EcophysiologyParms/GenusParmsQuery.aspx") 
      # which element
      opt1 <- remDr$findElement(using='xpath', value= '//*[@id="DropDownList3"]')
      # click
      opt1$clickElement()
      # we place a Family between the ''...
      opt2 <- remDr$findElement(using = 'xpath', paste(c(substr(element, 1, 21), 
                                                         substr(element, 22,23)), 
                                                       collapse = taxa.s[taxon]))
      # ...then click
      opt2$clickElement()
      # we select the Process Query button and click
      opt3 <- remDr$findElement(using = 'xpath', value= '//*[@id="Button1"]')
      opt3$clickElement()
      
      # ta-dan! we are inside the Genus Parameters page
      
      # I do not know why this works but it gives a list with all tables
      find.tables <- remDr$findElement(using = 'xpath', value= '//*[@id="DetailsView2"]')
      tables = find.tables$getPageSource()[[1]] %>% read_html() %>% html_table(fill = TRUE)
      
      # we select the tables with taxon parameters and combine them
      genus.list[[taxon]] = rbind(tables[[3]],tables[[4]],tables[[5]]) 
      
      # name columns after 1st row ("Genus","(Genus name)")
      genus.list[[taxon]] = genus.list[[taxon]] %>% 
        setNames(genus.list[[taxon]][1,]) %>% 
        # then drop that row
        dplyr::slice(-1)
      
      cat('\014')
      cat(paste0(round((taxon/length(taxa.s))*100), '%'))
      #Sys.sleep(.05)
      if (taxon == length(taxa.s)) cat('- Genera query complete
                                       ')
      
    }
    
    taxa.gen = genus.list %>% reduce(full_join, by = "Genus")
    taxa.gen = as.data.frame(taxa.gen)
    rownames(taxa.gen) = taxa.gen[,1]
    taxa.gen = taxa.gen[,-1, drop=FALSE]
    all.gen = as.data.frame(t(taxa.gen))
    all.gen[,3:26] = data.frame(lapply(all.gen[,3:26], as.numeric))
    all.gen = dplyr::rename(all.gen, N=Ngenus)
    names(all.gen) = gsub("Genavg", "Avg", names(all.gen))
  }else{
    # in case your input contains no genera
    all.gen = vector(mode = "list", length(taxa.s))
    print("No genera in input")
  }
  
  if (complete == TRUE) {
    taxa.all = dplyr::bind_rows(all.fam, all.gen)
  }else{
    taxa.all = dplyr::bind_rows(all.fam, all.gen)
    taxa.all = taxa.all[c("cp_value","feeding","N","AvgMass","StderrMass")]
  }
  
  
  # closes firefox remote browser (which you don't see if headless)
  remDr$closeWindow()
  
  # kills java
  system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)

  
  
  print("These taxa were not found in the database")
  print(setdiff(taxa, c(g.present,f.present)))
  return(taxa.all)
  
}



# EI = 100 * [e / (e+b)] 
# e representing the enrichment component, calculated as the weighted frequencies
# of Ba1 and Fu2 nematodes and b representing the basal food web component, 
# calculated as the weighted frequencies of Ba2 and Fu2 nematodes 


# Using, as a reference, the fecundity and life-course
# characteristics of the Ba2 guild, with a weighting of
# 0.8 (Fig. 1), the Ba1 guild is more fecund and has a
# shorter life course. Under available resource conditions, 
# the population increase rate is about four times
# as great as that of the Ba2 guild (Ferris et al., 1996a,b).
# Consequently, they are assigned a weighting of 3.2

# b=(Ba2+Fu2)*W2, where W2.= 0.8,
# e=(Ba1*W1)+(Fu2*W2), where W1 = 3.2 and W2.= 0.8
# s=(Ban*Wn+Can*Wn+Fun*Wn+Omn*Wn)
# where n=3-5, W3.= 1.8, W4.= 3.2, W5.= 5.0

Enrichment <- function(df,nemaplex) {
  #  e = 3.2*Ba1 + .8*Fu2
  #  b = .8*(Ba2 + Fu2)
  #  EI = 100*(e/(e+b))
  
  # I hate that this unintuitive mostrosity works like a charm
  # in english: "select those columns in df, whose names match the names of those
  # rows in nemaplex, where "cp_value" is 1 and "feeding" is 3"
  Ba1 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 1 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  
  Ba2 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  
  Fu2 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
  
  index = df %>% mutate(EI = 100*((3.2*rowSums(.[Ba1]) + 
                                    .8*rowSums(.[Fu2])) /
                                  (3.2*rowSums(.[Ba1]) + 
                                    .8*rowSums(.[Fu2]) + .8*(rowSums(.[Ba2]) + 
                                                             rowSums(.[Fu2])))),
                        .keep = "none")
  
  tryme2 = df %>% select_if(~!is.numeric(.x))
  
  out = cbind(tryme2,index)
  return(out)
}


Structure <- function(df, nemaplex) {  # This needs to be checked against Marcels values!!!
  
  # the guilds indicative of structure are:
  Ba3 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  Ba4 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  Ba5 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  
  Fu3 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
  Fu4 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
  Fu5 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
  
  Om3 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 8)
  Om4 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 8)
  Om5 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 8)
  
  Ca2 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5)
  Ca3 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 3 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5)
  Ca4 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 4 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5)
  Ca5 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 5 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 5)
  
  # the guilds indicative of basal characteristics are:
  Ba2 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  
  Fu2 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)
  
  # calculate the index
  index = df %>% mutate(SI = 100*( .8* rowSums(.[Ca2]) + # in the book this is not in
                                  1.8*(rowSums(.[Ba3]) + 
                                       rowSums(.[Fu3]) + 
                                       rowSums(.[Om3]) + 
                                       rowSums(.[Ca3]))+
                                  3.2*(rowSums(.[Ba4]) + 
                                       rowSums(.[Fu4]) + 
                                       rowSums(.[Om4]) + 
                                       rowSums(.[Ca4]))+
                                    5*(rowSums(.[Ba5]) + 
                                       rowSums(.[Fu5]) + 
                                       rowSums(.[Om5]) + 
                                       rowSums(.[Ca5]))) / 
                                 ( .8* rowSums(.[Ca2]) + # again
                                  1.8*(rowSums(.[Ba3]) + 
                                       rowSums(.[Fu3]) + 
                                       rowSums(.[Om3]) + 
                                       rowSums(.[Ca3]))+
                                  3.2*(rowSums(.[Ba4]) + 
                                       rowSums(.[Fu4]) + 
                                       rowSums(.[Om4]) + 
                                       rowSums(.[Ca4]))+
                                    5*(rowSums(.[Ba5]) + 
                                       rowSums(.[Fu5]) + 
                                       rowSums(.[Om5]) + 
                                       rowSums(.[Ca5])) + .8*(rowSums(.[Ba2]) + 
                                                              rowSums(.[Fu2]))), 
                        .keep = "none")
  
  tryme2 = df %>% select_if(~!is.numeric(.x))
  
  out = cbind(tryme2, index)
  return(out)
}

Channel <- function(df, nemaplex) { # This needs to be checked against Marcels values!!!
  # channel index (CI), 100 * (0.8Fu2/(3.2Ba1 + 0.8Fu2))
  
  Ba1 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 1 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 3)
  
  Fu2 = which(nemaplex$cp_value[match(colnames(df), rownames(nemaplex))] == 2 & 
                nemaplex$feeding[match(colnames(df), rownames(nemaplex))] == 2)

  index = df %>% mutate(CI = 100*(.8*rowSums(.[Fu2])/(3.2*rowSums(.[Ba1]) + 
                                                       .8*rowSums(.[Fu2]))),
                        .keep = "none")
  
  tryme2 = df %>% select_if(~!is.numeric(.x))
  
  out = cbind(tryme2, index)
  return(out)
}


all.indices <- function(df, nemaplex) {
  out = cbind(Enrichment(df, nemaplex),
               Structure(df, nemaplex),
                 Channel(df, nemaplex))
  return(out)
}