# Focused checks for the reviewer-related 0.1.1 update; base R only.
quiet <- function(fun, ...) {
  capture.output(value <- suppressWarnings(fun(...)))
  value
}
k <- 5
plan <- quiet(SampleSurvival::SampleSurvival, k=k, EPV=20, event_rate=0.18)
stopifnot(plan$n_total == 556, plan$events_needed == 100,
          plan$n_parameters == 5, plan$n_variables == plan$n_parameters)
verification <- quiet(SampleSurvival::VerifyEPV, n_final=plan$n_total,
                      n_events=plan$events_needed, k=k)
stopifnot(verification$EPV==20, verification$k==5)
stopifnot(inherits(try(quiet(SampleSurvival::VerifyEPV,n_final=300,
                           n_events=100,k=0),silent=TRUE),"try-error"))
stopifnot(inherits(try(quiet(SampleSurvival::VerifyEPV,n_final=300,
                           n_events=100,k=2.5),silent=TRUE),"try-error"))
for (lang in c("en","es")) {
  for (epv in c(5,10,20,30,50)) {
    output <- capture.output(ans <- SampleSurvival::VerifyEPV(
      n_final=500,n_events=epv*k,k=k,language=lang))
    stopifnot(ans$EPV==epv,
      !any(grepl("Excellent|Excelente|Acceptable|Aceptable|very unstable|muy inestables|has robust|tiene estimaciones robustas",output)))
    stopifnot(any(grepl(if(lang=="en") "excluding intercept" else "sin intercepto", output)),
              any(grepl(if(lang=="en") "does not assess" else "no evalua", output)))
  }
}
example_file <- system.file("examples","revision_examples.R",package="SampleSurvival")
stopifnot(nzchar(example_file))
invisible(capture.output(source(example_file,local=new.env())))
cat("Reviewer update checks passed for SampleSurvival 0.1.1\n")
