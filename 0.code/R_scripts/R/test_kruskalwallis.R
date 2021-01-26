test_kruskalwallis <- function(input){
  
  #shape dataframe
  df <- lapply(input, melt)
  df <- ldply(df, data.frame)
  
  #With brain regions as variables
  #ignoring stages - just overall means difference
  test <- kruskal.test(value ~ variable,
                       data = df)
  print(test)
  
  #post <-  dunnTest(value ~ variable,
  #                  data=df,
  #                  method="bh")   
  
  #In case you want to include it, requires FSA library
  #print(post)
  

}