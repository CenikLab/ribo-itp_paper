

make_name_mapper = function(experiment_names){
  exp_names        = sort( experiment_names )
  unique_raw_names = unique(exp_names )
  group_names      = sort(unique(unlist( lapply(strsplit( unique_raw_names, split = "-" ) , "[[", 3 ) ) ) )
  
  original_exp_names     = c()
  new_experiment_names   = c()
  list_counter           = 1
  experiment_name_mapper = list()
  
  for(g in group_names){
    elements_of_group = unique_raw_names[grep(  g, unique_raw_names )]
    size_of_group     = length( elements_of_group  )
    
    for( i in seq(1:size_of_group) ){
      this_basename = strsplit( elements_of_group[i], split = "-"   )[[1]][3]
      new_name = paste(  this_basename, i , sep = "-" )
      print( c( elements_of_group[i], new_name )  )
      
      original_exp_names   = c(original_exp_names, elements_of_group[i])
      new_experiment_names = c(new_experiment_names, new_name) 
      
      experiment_name_mapper[[list_counter]] = new_name
      list_counter = list_counter + 1
    }
    
  }
  
  names(  experiment_name_mapper ) = original_exp_names
  
  return(experiment_name_mapper)
}



