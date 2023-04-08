#' Connect to Oracle
#'
#' Connect to set up a connection to Oracle using RODBC and getPass
#'
#' @param schema Connection. Default = 'AFSC'
#' @export

get_connected <- function(schema='AFSC'){(echo=FALSE)
  username <- getPass::getPass(msg = "Enter your ORACLE Username: ")
  password <- getPass::getPass(msg = "Enter your ORACLE Password: ")
  channel  <- RODBC::odbcConnect(paste(schema),paste(username),paste(password), believeNRows=FALSE)
}