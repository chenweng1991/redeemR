
#' add_changes
#'
#' @param variant given a variant, output the changes 
#' @return changes
#' @export
add_changes <-function(variant){
    changes<-sub("^\\d+_", "",variant)
    return(changes)
}

#' add_types
#'
#' @param changes given a changes, output the type (transition or transversion)
#' @return changes
#' @export
add_types <- function(changes){
    types<- ifelse(changes %in% c("C_T","G_A","T_C","A_G"),"transition","transversion")
    return(types)
}

#' Annotate_base_change
#' 
#' @param  input redeem object,  it takes V.fitered, add the nucleotide change, 
#' @return redeem object with the V.fitered slot modified
Annotate_base_change <- function(redeem){
    redeem@V.fitered<- redeem@V.fitered %>% mutate(changes=add_changes(Variants)) %>% mutate(types=add_types(changes))
    return(redeem)
}
