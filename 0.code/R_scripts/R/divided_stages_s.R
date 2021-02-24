divided_stages_s <- function(permrun_s, returned_df){
  
  #divides permutation runs into stages
  fetal1 <- permrun_s[grep("153|150|113|103|149|114", permrun_s$variable), "variable"]
  fetal2 <- permrun_s[grep("178|154|B96|B97", permrun_s$variable), "variable"]
  fetal3 <- permrun_s[grep("B98|107|B92|159", permrun_s$variable), "variable"]
  birth_inf <- permrun_s[grep("155|194|121|132|139", permrun_s$variable), "variable"]
  inf_child <- permrun_s[grep("131|171|122|143|173", permrun_s$variable), "variable"]
  child <- permrun_s[grep("172|118|141|174|175", permrun_s$variable), "variable"]
  adolescence <- permrun_s[grep("124|119|105|127", permrun_s$variable), "variable"]
  adult <- permrun_s[grep("130|136|126|145|123|135", permrun_s$variable), "variable"]

}