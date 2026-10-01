s <- paste(readLines('tutorial_v2.qmd'),collapse='\n')
specs <- list(c('dat_moura2021','effect.size.id','vi','VCV',2),c('dat_lim2014','id','vi','vcv',2),c('dat_spain','effect_id','var_Hedges','VCV_spain',1))
n <- 0L
for(z in specs) {
 start <- paste0(z[1],'$',z[2],' <- factor(',z[1],'$',z[2],')')
 end <- paste0('as.numeric(',z[1],'$',z[3],')))')
 positions <- gregexpr(start,s,fixed=TRUE)[[1]]
 codes <- lapply(positions,function(i) {x<-substring(s,i);j<-regexpr(end,x,fixed=TRUE)[[1]];if(j<0||j>1500)return(NULL);substring(x,1,j+nchar(end)-1)})
 codes <- Filter(Negate(is.null),codes)
 stopifnot(length(codes)==as.integer(z[5]))
 for(code in codes) {
  for(ids in list(c('2','10','1'),factor(c('2','10','1'),levels=c('1','2','10')))) {
   e <- new.env();d<-data.frame(dummy=1:3);d[[z[2]]]<-ids;d[[z[3]]]<-c(.2,.1,.3);assign(z[1],d,e)
   eval(parse(text=code),e);m<-get(z[4],e)
   stopifnot(identical(as.numeric(diag(m)[as.character(ids)]),c(.2,.1,.3)))
  }
  e<-new.env();d[[z[2]]]<-c('1','1','2');assign(z[1],d,e)
  stopifnot(inherits(try(eval(parse(text=code),e),silent=TRUE),'try-error'))
  n <- n+1L
 }
}
stopifnot(n==5L)
cat('FIVE_BLOCKS_ALIGNMENT_PERMUTATION_AND_DUPLICATE_CHECKS_PASSED\n')
