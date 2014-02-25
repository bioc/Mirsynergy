test_mirsynergy <- function() {
  
  load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))
  
  # run mirsynergy clustering
  V2 <- mirsynergy(W, H, verbose=TRUE)

  # check with presaved module cluster assignment
  checkTrue(identical(V, V2))  
}