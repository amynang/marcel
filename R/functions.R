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


query_nemaplex = function(taxa) {
  require(tidyverse)
  require(rvest)
  require(rJava)
  require(RSelenium)
  
  # start the thingy
  driver <- rsDriver(browser = "firefox")
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
    
    
  }
  
  taxa.fam = family.list %>% reduce(full_join, by = "Family")
  taxa.fam = as.data.frame(taxa.fam)
  rownames(taxa.fam) = taxa.fam[,1]
  taxa.fam = taxa.fam[,-1]
  all.fam = as.data.frame(t(taxa.fam))
  
  
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
    
    
  }
  
  taxa.gen = genus.list %>% reduce(full_join, by = "Genus")
  taxa.gen = as.data.frame(taxa.gen)
  rownames(taxa.gen) = taxa.gen[,1]
  taxa.gen = taxa.gen[,-1]
  all.gen = as.data.frame(t(taxa.gen))
  
  taxa.all = list(all.gen, all.fam)
  
  
  
  
  print("These taxa were not found in the database")
  print(setdiff(taxa, c(g.present,f.present)))
  return(taxa.all)
  
  system("taskkill /im java.exe /f", intern=FALSE, ignore.stdout=FALSE)
  
}




Enrichment <- function() {

}

Structure <- function() {

}

Channel <- function() {

}
