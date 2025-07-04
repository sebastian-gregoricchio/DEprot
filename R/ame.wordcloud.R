### motifs AME

ame.wordcloud =
  function(ame.results,
           e.value.column.id = "E-value",
           pattern.motifs.keep = "[.]0[.]",
           pattern.remove.from.ids = "[.].*",
           display.top.n = NA,
           remove.mouse = TRUE,
           palette = viridis::mako(n = 10, direction = -1, begin = 0.25),
           eval.scale.limits = c(NA,NA),
           font.range = c(2,15),
           font.breaks = NULL,
           wordcloud.shape = "circle",
           show.legend = TRUE,
           title = NULL) {
    
    ### Libraries
    require(ggplot2)
    require(ggwordcloud)
    require(dplyr)
    #require(ggtext)
    
    
    ### Read table
    eval.column = "character"
    names(eval.column) = e.value.column.id
    
    ame = read.delim(file = ame.results, comment.char = "#", colClasses = eval.column, header = T, sep = "\t", check.names = F)
    colnames(ame)[which(colnames(ame) == e.value.column.id)] = "E.value"
    
    
    ### Filter table and compute -log10(Eval)
    ## too big exp will be converted to zeros. So, Eval is imported as character,
    ## split in num and exp.base, and then logs computed separately -log10(Xe-Y) = Y - log10(X)
    
    ame = 
      ame %>%
      dplyr::filter(grepl(pattern.motifs.keep, motif_ID)) %>%
      dplyr::mutate(motif_ID = gsub(pattern.remove.from.ids,"",motif_ID)) %>%
      tidyr::separate(col = E.value, into = c("e.num", "e.exp"), sep = "e-", remove = F, convert = T, fill = "right") %>%
      dplyr::mutate(e.exp = ifelse(is.na(e.exp), yes = 0, no = e.exp)) %>%
      dplyr::mutate(m.log10.eval = e.exp - log10(e.num)) %>%
      dplyr::arrange(as.numeric(rank))
    
      
    ## remove mouse motifs (if required)
    if (remove.mouse == TRUE) {ame = ame %>% dplyr::filter(toupper(motif_ID) == motif_ID)}
    
    
    ## Define how many motifs to display
    if (!is.na(display.top.n)) {
      if (display.top.n >= nrow(ame)) {
        n.motifs = nrow(ame)
      } else {
        n.motifs = display.top.n
      }
    } else {
      n.motifs = nrow(ame)
    }
    
    ame.plot = ame[1:n.motifs,]
    
    
    
    ### generate wordcloud
    wordcloud =
      ggplot(data = ame.plot,
             aes(label = motif_ID,
                 size = m.log10.eval,
                 color = m.log10.eval)) +
      ggwordcloud::geom_text_wordcloud(eccentricity = 1,
                                       shape = wordcloud.shape,
                                       show.legend = show.legend,
                                       rm_outside = T) +
      scale_color_gradientn(name = "-log<sub>10</sub>(*E-value*)",
                            colours = palette,
                            limits = eval.scale.limits) +
      ggtitle(title) +
      theme(panel.background = element_blank(),
            legend.text = ggtext::element_markdown(),
            legend.title = ggtext::element_markdown(),
            plot.title = ggtext::element_markdown(hjust = 0.5))
    
    if (is.null(font.breaks)) {
      wordcloud = 
        wordcloud +
        scale_size(name = "-log<sub>10</sub>(*E-value*)",
                   limits = eval.scale.limits,
                   range = font.range)
    } else {
      wordcloud = 
        wordcloud +
        scale_size(name = "-log<sub>10</sub>(*E-value*)",
                   limits = eval.scale.limits,
                   range = font.range,
                   breaks = font.breaks)
    }
    
    
    # return plot
    return(wordcloud)
    
  } # END function