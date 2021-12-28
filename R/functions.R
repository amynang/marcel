
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

Enrichment <- function() {

}

Structure <- function() {

}

Channel <- function() {

}
