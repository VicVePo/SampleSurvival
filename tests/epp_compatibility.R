# EPP interface and legacy compatibility checks; base R only.
quiet <- function(fun, ...) {
  capture.output(value <- suppressWarnings(fun(...)))
  value
}
calc <- SampleSurvival::SampleSurvival
verify <- SampleSurvival::VerifyEPP
legacy_verify <- SampleSurvival::VerifyEPV
stopifnot(identical(formals(verify),formals(legacy_verify)))
for (lang in c('en','es')) {
  for (target in c(5,10,20,30,50)) {
    old <- quiet(calc,k=5,event_rate=0.18,EPV=target,language=lang)
    text <- capture.output(new <- calc(k=5,event_rate=0.18,EPP=target,language=lang))
    stopifnot(identical(old,new),new$target_EPP==target,new$target_EPV==target,
              any(grepl('EPP',text)),!any(grepl('EPV',text)))
    both <- quiet(calc,k=5,event_rate=0.18,EPV=target,EPP=target,language=lang)
    stopifnot(identical(new,both))
    a <- quiet(verify,n_final=500,n_events=target*5,k=5,language=lang)
    b <- quiet(legacy_verify,n_final=500,n_events=target*5,k=5,language=lang)
    stopifnot(identical(a,b),a$EPP==target,a$EPV==target)
  }
  text <- capture.output(sc <- calc(k=5,event_rate=0.18,scenarios=TRUE,language=lang))
  stopifnot(identical(sc$EPP,sc$EPV),any(grepl('EPP',text)),!any(grepl('EPV',text)))
}
stopifnot(quiet(calc,k=5,event_rate=0.18)$target_EPP==20)
stopifnot(inherits(try(quiet(calc,k=5,event_rate=0.18,EPV=10,EPP=20),silent=TRUE),'try-error'))
for (bad in list(0,-1,NA_real_,Inf,numeric(),c(10,20),'20')) {
  stopifnot(inherits(try(quiet(calc,k=5,event_rate=0.18,EPP=bad),silent=TRUE),'try-error'))
}
cat('EPP interface and compatibility checks passed for SampleSurvival 0.1.2\n')
