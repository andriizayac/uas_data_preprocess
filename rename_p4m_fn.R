rename_p4m <- function(f) {
	lapply(f, function(x){
	s <- str_split(x, pattern = " ")
	if(lengths(s) == 1){
		x_new <- gsub("_", "_1", s) 
	} else {
		x_new <- gsub("_", paste0("_", str_sub(s[[1]][2], 2,2)), paste0(s[[1]][1], ".TIF"))
	}
	return(x_new)
	}) |> unlist() -> out
	
	return(out)
}