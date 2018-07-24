#'@title GDM variable selection
#'
#'@description Perform backward elimination to determine variable significance in a gdm 
#'
#'@param data (string, data.frame or matrix). Input table, can be in memory (data.frame or matrix) or filepath to .CSV
#'@param dst (string) Filepath to write outputs to. A subdir (taking the name of outname below) will be created here.
#'@param outname (string) Filename to use for writing outputs. Default: 'gdm_variable_selection'
#'@param permutations (int) Number of permutations to run for each variable. Default: 100.
#'@param geo (boolean) Whether or not to include geographic distance as a predictor. Default: TRUE
#'@param load_output (boolean) If TRUE (default) output files will be loaded into memory in a list.
#'@param verbose (boolean) Print out messages to console. Default: FALSE
#'
#'@return Files written to dst.
#'
#'@examples gdm.variable_select(inputTable, dst_dir, analysis_1, 50, TRUE, FALSE, TRUE)
#'
#'@useDynLib gdmToolbox SaveGDMParams ExtractAndUpdateQuantilesSigTest DoSigTestGDM
#'
#'@importFrom gdm gdm
#'@importFrom assertthat assert_that is.string
#'@importFrom data.table fread
#'
#'@export

gdm.variable_select = function(data, dst, outname = 'gdm_variable_selection',
                               permutations = 100, geo = TRUE, load_output = TRUE,
                               verbose = FALSE){
  
  ## check what data is
  data_class = class(data)
  in_memory = sum(unlist(lapply(c('data.frame', 'matrix', 'data.table'), 
                         function(x) grep(x, data_class))))
  if(in_memory > 0){
    ## in memory is true
    data = data.frame(data)
    numpreds = (ncol(data) - 6) / 2
    ## write to tmp
    tmp_loc = tempfile(fileext = '.csv')
    write.csv(data, tmp_loc, row.names = FALSE)
    data = tmp_loc
  
  } else {
    ## probably filepath - double check
    if(data_class == 'character'){
      ## assume filepath
      if (!file.exists(data)){
        stop('data is of class character but is not a valid file path')
      } else {
        numpreds = (ncol(fread(data, nrows = 1L)) - 6) / 2
      }
    } else {
      stop('Cannot determine class of input data')
    }  
  }
  
  ## dst arg must exist
  assert_that(is.notempty.string(dst))
  
  ## check dst is not empty
  try_dst = paste0(dst, '/', outname)
  if(file.exists(try_dst)){
    warning(cat('Destination folder ', outname, ' already exists -\n', 
                  'a new sub folder will be created in this folder',
            sep = ''))
    dst = try_dst
  } else {
    dst = try_dst
    dir.create(dst)
  }
  
  ## check permutations arg
  assert_that(is.numeric(permutations))
  
  ## geo
  assert_that(is.logical(geo))  

  if (verbose){
    if (geo){
      print("Including geographic distance as a predictor")
    } else {
      print("Excluding geographic distance as a predictor")
    } 
  }
  
  ## find dll
  this_lib = dirname(find.package('gdmToolbox', .libPaths()))
  ## check architecture and then guess filepath to dll - possibly not 
  ## a solid way of doing this...
  arch = Sys.info()[['release']]
  opts = c('86', '64')
  check_arch = lapply(opts, function(x) grep(x, arch))
  check_arch = unlist(lapply(lapply(check_arch, length), function(x) x == 1))
  arch = opts[check_arch]
  if(arch == '64'){
    dllpath = paste0(this_lib, '/gdmToolbox/libs/x64/gdmToolbox.dll')
  }
  if(arch == '86'){
    dllpath = paste0(this_lib, '/gdmToolbox/libs/i386/gdmToolbox.dll')
  }
  assert_that(file.exists(dllpath))
  ## format path and load
  dllpath = paste(strsplit(dllpath, '/')[[1]], collapse = '\\')
  dyn.load(dllpath)
  
  ## dll args to write param file
  wdpath = gsub( "\\\\",  "/", dst)
  wdpath = gsub('/', '\\\\', wdpath)
  paramFilePath = paste0(outname, '.txt')
  datatable = gsub( "\\\\",  "/", data)
  datatable = gsub('/', '\\\\', data)
  numpreds = as.numeric(numpreds)
  do_geo = geo
  
  z0 <- .C( "SaveGDMParams", 
            wdpath, 
            paramFilePath, 
            datatable, 
            as.integer(numpreds), 
            as.integer(do_geo), PACKAGE = 'gdmToolbox')
  
  fullparamfilepath <- paste0(wdpath, "\\", paramFilePath)
  z1 <- .C( "ExtractAndUpdateQuantilesSigTest", fullparamfilepath, 
            PACKAGE = 'gdmToolbox')
  
  z2 <- .C( "DoSigTestGDM", fullparamfilepath, as.integer(permutations), 
            PACKAGE = 'gdmToolbox')
  
  dyn.unload(dllpath)
  
  cat('Significance testing done. Output files located here:\n', dst, sep = '')
  
  if(load_output){
    f = list.files(dst, pattern = '.csv')
    outputs = lapply(f, read.csv)
    names(outputs) = gsub('.csv', '', unlist(lapply(f, basename)))
    return(outputs)
  }
              
  
}



#'@title Helper
#'
#'@export

is.notempty.string = function(x) {
  is.string(x) && !is.na(x) && nchar(x)>0
}

