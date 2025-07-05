#' Check redeem object consistency
#'
#' @description
#' Verifies that the counts matrix and binary counts matrix in a redeem object have matching dimensions and identical row and column names.
#'
#' @param obj Redeem object with slots \code{Cts.Mtx} and \code{Cts.Mtx.bi}.
#'
#' @return Prints confirmation messages for each consistency check.
#' @examples
#' check_redeem(my_redeem)
#' @export
check_redeem <-function(obj){
    if (dim(obj@Cts.Mtx)[1]==dim(obj@Cts.Mtx.bi)[1] & dim(obj@Cts.Mtx)[2]==dim(obj@Cts.Mtx)[2]){
        print("confirmed: Mtx and Mtx.bi dimensions are the same")
    }
    else{
        print("dimensions are different, double check")
    }
    if(all(row.names(obj@Cts.Mtx)== row.names(obj@Cts.Mtx.bi))){
        print("confirmed: all cell names match")
    }
    if(all(colnames(obj@Cts.Mtx)== colnames(obj@Cts.Mtx.bi))){
        print("confirmed: all variant names match")
    }
    # print("transversion and transition")
    # obj@V.fitered$types %>% table %>% prop.table
}