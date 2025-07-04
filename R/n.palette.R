# generate palette gradient
n.palette = 
  function(palette.name,
           n.colors = 101) {
    
    ## define palette sets
    pal.9 = c("YlOrRd","YlOrBr","YlGnBu","YlGn","Reds","RdPu","Purples","PuRd","PuBuGn","PuBu",
              "OrRd","Oranges","Greys","Greens","GnBu","BuPu","BuGn","Blues")
    pal.11 = c("Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PuOr", "PRGn", "PiYG")
    
    ## check palette name
    if (!(tolower(palette.name) %in% tolower(c(pal.9,pal.11))) | !(class(palette.name) %in% "character")) {
      return(warning(paste0("The plette chosen is not available, please chose one among: ",
                            "'",paste0(pal.9, collapse = "', '"),
                            "'",paste0(pal.11, collapse = "', '"),"'.")))
    } else {
      # get the name with the correct capitalization
      name = c(pal.9,pal.11)[tolower(palette.name) == tolower(c(pal.9,pal.11))]
      
      ## Define the palette size
      palette.size = ifelse(test = name %in% pal.9, yes = 9, no = 11)
    }
    
    return(colorRampPalette(RColorBrewer::brewer.pal(n = palette.size, name = name))(n.colors))
  }



