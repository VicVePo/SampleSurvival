# Version 0.1.1: predictor-parameter counting example.
# This grid defines coding; it is not a patient dataset.
design <- expand.grid(age=c(30,40,50,60,70),
                      sex=factor(c("F","M")),
                      education=factor(c("A","B","C","D")))
mm <- stats::model.matrix(~ age + sex + education, data=design)
k <- sum(colnames(mm) != "(Intercept)")
stopifnot(k == 5L, qr(mm)$rank == ncol(mm))
planned <- SampleSurvival::SampleSurvival(k=k, EPV=20, event_rate=0.18)
SampleSurvival::VerifyEPV(n_final=planned$n_total,
                       n_events=planned$events_needed, k=k)

# A log transformation uses one coefficient; nonlinear terms can use more.
forms <- list(linear=~age+sex+education,
              logarithmic=~log(age)+sex+education,
              quadratic=~age+I(age^2)+sex+education,
              spline=~splines::ns(age,df=3)+sex+education,
              interaction=~age*sex+education)
parameter_counts <- vapply(forms, function(f) {
  m <- stats::model.matrix(f, data=design)
  sum(colnames(m) != "(Intercept)")
}, integer(1))
stopifnot(identical(unname(parameter_counts),c(5L,5L,6L,7L,6L)))
print(parameter_counts)

# Reuse the computed target rather than typing 1333 or 1334 manually.
renal <- SampleSurvival::SampleSurvival(k=12, event_rate=0.18, EPV=20)
logistics <- SampleSurvival::SurvivalLogistics(
    n_final=renal$n_total, loss_rate=0.12, eligibility_rate=0.75,
    subjects_per_month=25, follow_up_months=60, annual_event_rate=0.036)
stopifnot(renal$n_total==1334, logistics$initial_cohort==1516,
          logistics$people_to_invite==2022, logistics$total_duration_months==141)
