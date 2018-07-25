#'@title GDM variable selection
#'
#'@description Provides an automated way to select a single model using gdm.variable_significance output
#'
#'@param data (string, data.frame or matrix). Input table, can be in memory (data.frame or matrix) or filepath to .CSV
#'@param dst (string) Filepath to write outputs to. A subdir (taking the name of outname below) will be created here. Default (NULL) will write to a temporary dir.
#'@param outname (string) Filename to use for writing outputs. Default (NULL) will assign name of 'gdm_variable_significance.'
#'@param ... Other args to be passed to \code{\link{gdm.variable_significance}}
#'
#'@export

gdm.variable_selection = function(data, dst = NULL, outname = NULL, ...){
  if (is.null(dst)) dst = tempdir()
  if (is.null(outname)) outname = 'gdm_variable_significnace'
  
  pass_args = list(
    data = data, 
    dst = dst,
    outname = outname,
    load_output = TRUE
  )
  
  user_args = list(...)
  if (length(user_args)) {
    opts = formals(gdm.variable_significance)
    for (i in 1:length(user_args)) {
      this_opt <- names(user_args)[i]
      if (! (this_opt %in% names(opts))) {
        cat(paste0("'", this_opt, "'", ' is not a valid option. Should be one of:'), 
            names(opts), sep = '\n')
      } else {
        pass_args[[this_opt]] = paste(user_args[i])
      }
    } # end user_opts 
  } 
  
  outputs = do.call(gdm.variable_significance, pass_args)
  
  ## do something with outputs...
  ## just return them for now
  return(outputs)
  
}
