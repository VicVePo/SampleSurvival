#' Sample Size Calculation for Survival Analysis (Cox Regression)
#' 
#' Calculates the required sample size for multivariable analysis of associated
#' factors in open cohort studies (survival analysis). Uses the events per 
#' variable (EPV) method appropriate for Cox regression estimating hazard 
#' ratios (HR).
#' 
#' @param k Number of independent variables to include in the multivariable model
#' @param event_rate Expected proportion that will experience the event during 
#'   follow-up (between 0 and 1). For example, 0.15 = 15% will have the event
#' @param EPV Events per variable (recommended 20-50)
#' @param scenarios Logical. If TRUE, calculates sample sizes for multiple 
#'   EPV values (10, 20, 30, 40, 50)
#' @param language Language for messages: 'en' (English) or 'es' (Spanish). Default is 'en'
#' @return List or data.frame with results
#' @export
#' @examples
#' # Open cohort with 20% expected events
#' SampleSurvival(k = 10, event_rate = 0.20, EPV = 20)
#' 
#' # Spanish version
#' SampleSurvival(k = 10, event_rate = 0.20, EPV = 20, language = 'es')
#' 
#' # View multiple scenarios
#' SampleSurvival(k = 10, event_rate = 0.20, scenarios = TRUE)
SampleSurvival <- function(k, event_rate, EPV = 20, scenarios = FALSE, language = 'en') {
  
  if (!language %in% c('en', 'es')) {
    stop('language must be "en" or "es"')
  }
  
  if (missing(k)) {
    stop(ifelse(language == 'es',
                'Debe especificar k (numero de variables)',
                'You must specify k (number of variables)'))
  }
  if (missing(event_rate)) {
    stop(ifelse(language == 'es',
                'Debe especificar event_rate esperada',
                'You must specify expected event_rate'))
  }
  if (k <= 0 || k != round(k)) {
    stop(ifelse(language == 'es',
                'k debe ser entero positivo',
                'k must be a positive integer'))
  }
  if (event_rate <= 0 || event_rate >= 1) {
    stop(ifelse(language == 'es',
                'event_rate debe estar entre 0 y 1 (ejemplo: 0.15 para 15%)',
                'event_rate must be between 0 and 1 (example: 0.15 for 15%)'))
  }
  if (EPV <= 0) {
    stop(ifelse(language == 'es',
                'EPV debe ser positivo',
                'EPV must be positive'))
  }
  
  if (scenarios) {
    valores_EPV <- c(10, 20, 30, 40, 50)
    res <- data.frame(
      EPV = valores_EPV,
      eventos_necesarios = valores_EPV * k,
      censurados_esperados = ceiling(((valores_EPV * k) / event_rate) - (valores_EPV * k)),
      n_total = ceiling((valores_EPV * k) / event_rate)
    )
    
    if (language == 'es') {
      names(res) <- c('EPV', 'eventos_necesarios', 'censurados_esperados', 'n_total')
      cat('\nESCENARIOS DE TAMAÑO MUESTRAL - ANALISIS SUPERVIVENCIA\n')
      cat('Variables (k):', k, '\n')
      cat('Tasa de eventos esperada:', event_rate*100, '%\n\n')
      print(res, row.names = FALSE)
      cat('\nRecomendacion: EPV >= 20\n')
      cat('Nota: n_total incluye tanto eventos como censurados\n\n')
    } else {
      names(res) <- c('EPV', 'events_needed', 'censored_expected', 'n_total')
      cat('\nSAMPLE SIZE SCENARIOS - SURVIVAL ANALYSIS\n')
      cat('Variables (k):', k, '\n')
      cat('Expected event rate:', event_rate*100, '%\n\n')
      print(res, row.names = FALSE)
      cat('\nRecommendation: EPV >= 20\n')
      cat('Note: n_total includes both events and censored\n\n')
    }
    return(invisible(res))
  }
  
  if (EPV < 10) {
    warning(ifelse(language == 'es',
                   'EPV < 10 muy bajo. Se recomienda EPV >= 20',
                   'EPV < 10 very low. EPV >= 20 is recommended'))
  }
  
  # Calculations
  events_needed <- k * EPV
  n_total <- ceiling(events_needed / event_rate)
  n_censored <- n_total - events_needed
  
  resultados <- list(
    design = ifelse(language == 'es', 
                    'Cohorte abierta (analisis de supervivencia)',
                    'Open cohort (survival analysis)'),
    model = ifelse(language == 'es', 'Regresion de Cox', 'Cox regression'),
    association_measure = 'Hazard Ratio (HR)',
    n_variables = k,
    target_EPV = EPV,
    expected_event_rate = event_rate,
    events_needed = events_needed,
    censored_expected = n_censored,
    n_total = n_total,
    language = language
  )
  class(resultados) <- c('SampleSurvival', 'list')
  
  if (language == 'es') {
    cat('\n=== TAMAÑO MUESTRAL - ANALISIS SUPERVIVENCIA (COX) ===\n')
    cat('Modelo: Regresion de Cox (HR)\n')
    cat('Variables (k):', k, '\n')
    cat('EPV:', EPV, '\n')
    cat('Tasa de eventos esperada:', event_rate*100, '%\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Eventos necesarios:', events_needed, '\n')
    cat('Censurados esperados:', n_censored, '\n')
    cat('\n>>> TAMAÑO TOTAL:', n_total, '<<<\n\n')
    
    if (EPV < 20) {
      cat('⚠️  ADVERTENCIA: EPV < 20 puede comprometer validez\n')
      cat('   de estimaciones en modelos de Cox\n\n')
    } else if (EPV >= 50) {
      cat('✓ Excelente: EPV >= 50 proporciona estimaciones muy robustas\n\n')
    }
    
    if (event_rate < 0.10) {
      cat('ℹ️  Nota: Tasa de eventos baja (<10%). El seguimiento debe ser\n')
      cat('   lo suficientemente largo para alcanzar los eventos necesarios\n\n')
    } else if (event_rate > 0.50) {
      cat('ℹ️  Nota: Tasa de eventos alta (>50%). Considere si el supuesto\n')
      cat('   de riesgos proporcionales es apropiado\n\n')
    }
  } else {
    cat('\n=== SAMPLE SIZE - SURVIVAL ANALYSIS (COX) ===\n')
    cat('Model: Cox regression (HR)\n')
    cat('Variables (k):', k, '\n')
    cat('EPV:', EPV, '\n')
    cat('Expected event rate:', event_rate*100, '%\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Events needed:', events_needed, '\n')
    cat('Censored expected:', n_censored, '\n')
    cat('\n>>> TOTAL SIZE:', n_total, '<<<\n\n')
    
    if (EPV < 20) {
      cat('⚠️  WARNING: EPV < 20 may compromise validity\n')
      cat('   of estimates in Cox models\n\n')
    } else if (EPV >= 50) {
      cat('✓ Excellent: EPV >= 50 provides very robust estimates\n\n')
    }
    
    if (event_rate < 0.10) {
      cat('ℹ️  Note: Low event rate (<10%). Follow-up must be\n')
      cat('   long enough to reach the needed events\n\n')
    } else if (event_rate > 0.50) {
      cat('ℹ️  Note: High event rate (>50%). Consider if the\n')
      cat('   proportional hazards assumption is appropriate\n\n')
    }
  }
  
  return(invisible(resultados))
}

#' Logistical Planning for Survival Studies
#' 
#' Calculates logistical requirements for a survival study, considering
#' recruitment, follow-up, and losses.
#' 
#' @param n_final Required sample size
#' @param loss_rate Proportion of loss to follow-up (0-1)
#' @param eligibility_rate Proportion of eligible people (0-1)
#' @param subjects_per_month People that can be recruited per month
#' @param follow_up_months Maximum follow-up time in months
#' @param annual_event_rate Annual event rate (to estimate when needed events will be reached)
#' @param language Language for messages: 'en' (English) or 'es' (Spanish). Default is 'en'
#' @return List with logistical requirements
#' @export
SurvivalLogistics <- function(n_final,
                             loss_rate,
                             eligibility_rate,
                             subjects_per_month,
                             follow_up_months,
                             annual_event_rate,
                             language = 'en') {
  
  if (!language %in% c('en', 'es')) {
    stop('language must be "en" or "es"')
  }
  
  if (missing(n_final)) {
    stop(ifelse(language == 'es',
                'Especifique n_final',
                'Specify n_final'))
  }
  if (missing(loss_rate)) {
    stop(ifelse(language == 'es',
                'Especifique loss_rate',
                'Specify loss_rate'))
  }
  if (missing(eligibility_rate)) {
    stop(ifelse(language == 'es',
                'Especifique eligibility_rate',
                'Specify eligibility_rate'))
  }
  if (missing(subjects_per_month)) {
    stop(ifelse(language == 'es',
                'Especifique subjects_per_month',
                'Specify subjects_per_month'))
  }
  if (missing(follow_up_months)) {
    stop(ifelse(language == 'es',
                'Especifique follow_up_months',
                'Specify follow_up_months'))
  }
  if (missing(annual_event_rate)) {
    stop(ifelse(language == 'es',
                'Especifique annual_event_rate',
                'Specify annual_event_rate'))
  }
  
  if (loss_rate < 0 || loss_rate >= 1) {
    stop(ifelse(language == 'es',
                'loss_rate debe estar entre 0 y menor a 1',
                'loss_rate must be between 0 and less than 1'))
  }
  if (eligibility_rate <= 0 || eligibility_rate > 1) {
    stop(ifelse(language == 'es',
                'eligibility_rate debe estar entre 0 y 1',
                'eligibility_rate must be between 0 and 1'))
  }
  
  # Adjust for loss to follow-up
  n_recruit <- ceiling(n_final / (1 - loss_rate))
  
  # Adjust for eligibility
  n_invite <- ceiling(n_recruit / eligibility_rate)
  
  # Recruitment time
  months_recruitment <- ceiling(n_invite / subjects_per_month)
  
  # Total study time
  total_months <- months_recruitment + follow_up_months
  
  # Expected losses
  n_losses <- n_recruit - n_final
  
  resultados <- list(
    final_cohort = n_final,
    initial_cohort = n_recruit,
    people_to_invite = n_invite,
    expected_losses = n_losses,
    recruitment_months = months_recruitment,
    follow_up_months = follow_up_months,
    total_duration_months = total_months,
    annual_event_rate = annual_event_rate,
    parameters = list(
      loss_rate = loss_rate,
      eligibility_rate = eligibility_rate,
      subjects_per_month = subjects_per_month
    ),
    language = language
  )
  class(resultados) <- c('SurvivalLogistics', 'list')
  
  if (language == 'es') {
    cat('\n=== LOGISTICA - ANALISIS SUPERVIVENCIA ===\n\n')
    cat('PARAMETROS:\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Tasa de perdida esperada:', loss_rate*100, '%\n')
    cat('Tasa de elegibilidad:', eligibility_rate*100, '%\n')
    cat('Capacidad reclutamiento:', subjects_per_month, 'personas/mes\n')
    cat('Tiempo maximo seguimiento:', follow_up_months, 'meses\n')
    cat('Tasa de eventos:', annual_event_rate*100, '% anual\n\n')
    
    cat('REQUERIMIENTOS:\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Cohorte final requerida:', n_final, 'participantes\n')
    cat('Cohorte inicial a reclutar:', n_recruit, 'participantes\n')
    cat('  (ajustado por', loss_rate*100, '% perdida)\n')
    cat('Personas a invitar:', n_invite, 'personas\n')
    cat('  (ajustado por', eligibility_rate*100, '% elegibilidad)\n')
    cat('Perdidas esperadas:', n_losses, 'participantes\n\n')
    
    cat('TIEMPO ESTIMADO:\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Reclutamiento:', months_recruitment, 'meses\n')
    cat('Seguimiento maximo:', follow_up_months, 'meses\n')
    cat('\n>>> DURACION TOTAL:', total_months, 'meses (',
        round(total_months/12, 1), 'años) <<<\n')
    cat('─────────────────────────────────────────────────\n\n')
    
    if (total_months > 60) {
      cat('⚠️  ADVERTENCIA: Estudio > 5 años. Evaluar factibilidad\n\n')
    } else if (total_months > 36) {
      cat('ℹ️  Nota: Estudio > 3 años. Planifique recursos a largo plazo\n\n')
    }
    
    if (loss_rate > 0.20) {
      cat('⚠️  ADVERTENCIA: Tasa de perdida alta (>20%)\n')
      cat('   La censura informativa puede sesgar resultados\n')
      cat('   Implemente estrategias rigurosas de retencion\n\n')
    }
    
    if (annual_event_rate < 0.05) {
      cat('ℹ️  Nota: Tasa de eventos baja (<5% anual)\n')
      cat('   Puede requerir seguimiento mas prolongado\n\n')
    }
  } else {
    cat('\n=== LOGISTICS - SURVIVAL ANALYSIS ===\n\n')
    cat('PARAMETERS:\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Expected loss rate:', loss_rate*100, '%\n')
    cat('Eligibility rate:', eligibility_rate*100, '%\n')
    cat('Recruitment capacity:', subjects_per_month, 'people/month\n')
    cat('Maximum follow-up time:', follow_up_months, 'months\n')
    cat('Event rate:', annual_event_rate*100, '% annual\n\n')
    
    cat('REQUIREMENTS:\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Final cohort required:', n_final, 'participants\n')
    cat('Initial cohort to recruit:', n_recruit, 'participants\n')
    cat('  (adjusted for', loss_rate*100, '% loss)\n')
    cat('People to invite:', n_invite, 'people\n')
    cat('  (adjusted for', eligibility_rate*100, '% eligibility)\n')
    cat('Expected losses:', n_losses, 'participants\n\n')
    
    cat('ESTIMATED TIME:\n')
    cat('─────────────────────────────────────────────────\n')
    cat('Recruitment:', months_recruitment, 'months\n')
    cat('Maximum follow-up:', follow_up_months, 'months\n')
    cat('\n>>> TOTAL DURATION:', total_months, 'months (',
        round(total_months/12, 1), 'years) <<<\n')
    cat('─────────────────────────────────────────────────\n\n')
    
    if (total_months > 60) {
      cat('⚠️  WARNING: Study > 5 years. Evaluate feasibility\n\n')
    } else if (total_months > 36) {
      cat('ℹ️  Note: Study > 3 years. Plan for long-term resources\n\n')
    }
    
    if (loss_rate > 0.20) {
      cat('⚠️  WARNING: High loss rate (>20%)\n')
      cat('   Informative censoring may bias results\n')
      cat('   Implement rigorous retention strategies\n\n')
    }
    
    if (annual_event_rate < 0.05) {
      cat('ℹ️  Note: Low event rate (<5% annual)\n')
      cat('   May require longer follow-up\n\n')
    }
  }
  
  return(invisible(resultados))
}

#' Post-Study EPV Verification (Survival Analysis)
#' 
#' Verifies if the observed EPV at the end of survival analysis is adequate
#' 
#' @param n_final Final sample in analysis
#' @param n_events Number of events observed (not censored)
#' @param k Number of variables in the Cox model
#' @param language Language for messages: 'en' (English) or 'es' (Spanish). Default is 'en'
#' @return List with observed EPV and evaluation
#' @export
VerifyEPV <- function(n_final, n_events, k, language = 'en') {
  
  if (!language %in% c('en', 'es')) {
    stop('language must be "en" or "es"')
  }
  
  if (missing(n_final)) {
    stop(ifelse(language == 'es',
                'Especifique n_final',
                'Specify n_final'))
  }
  if (missing(n_events)) {
    stop(ifelse(language == 'es',
                'Especifique n_events',
                'Specify n_events'))
  }
  if (missing(k)) {
    stop(ifelse(language == 'es',
                'Especifique k (numero de variables)',
                'Specify k (number of variables)'))
  }
  
  EPV_observed <- n_events / k
  observed_event_rate <- n_events / n_final
  n_censored <- n_final - n_events
  
  if (language == 'es') {
    cat('\n=== VERIFICACIÓN EPV POST-ESTUDIO (COX) ===\n')
    cat('Muestra final analizada:', n_final, '\n')
    cat('Eventos observados:', n_events, '\n')
    cat('Censurados:', n_censored, '\n')
    cat('Variables en modelo Cox:', k, '\n')
    cat('Tasa de eventos observada:', round(observed_event_rate*100, 2), '%\n')
    cat('\n>>> EPV OBSERVADO:', round(EPV_observed, 2), '<<<\n\n')
    
    if (EPV_observed < 10) {
      cat('❌ CRÍTICO: EPV < 10 en modelo de Cox\n')
      cat('   Los estimadores del HR son muy inestables\n')
      cat('   Recomendacion: Reducir numero de variables o considerar\n')
      cat('   tecnicas de penalizacion (ridge Cox, lasso Cox)\n\n')
    } else if (EPV_observed < 20) {
      cat('⚠️  ADVERTENCIA: EPV < 20 en modelo de Cox\n')
      cat('   Interpretacion cautelosa requerida\n')
      cat('   Recomendacion: Reportar como limitacion y realizar\n')
      cat('   analisis de sensibilidad. Verificar supuesto de\n')
      cat('   riesgos proporcionales cuidadosamente\n\n')
    } else if (EPV_observed >= 20 && EPV_observed < 30) {
      cat('✓ Aceptable: EPV >= 20\n')
      cat('  El modelo cumple el minimo recomendado para Cox\n\n')
    } else {
      cat('✓✓ Excelente: EPV >= 30\n')
      cat('   El modelo de Cox tiene estimaciones robustas\n\n')
    }
    
    # Additional warning about censoring
    censoring_proportion <- n_censored / n_final
    if (censoring_proportion > 0.50) {
      cat('ℹ️  Nota: Alta proporcion de censura (>', round(censoring_proportion*100), '%)\n')
      cat('   Verifique si la censura es informativa\n\n')
    }
  } else {
    cat('\n=== POST-STUDY EPV VERIFICATION (COX) ===\n')
    cat('Final sample analyzed:', n_final, '\n')
    cat('Events observed:', n_events, '\n')
    cat('Censored:', n_censored, '\n')
    cat('Variables in Cox model:', k, '\n')
    cat('Observed event rate:', round(observed_event_rate*100, 2), '%\n')
    cat('\n>>> OBSERVED EPV:', round(EPV_observed, 2), '<<<\n\n')
    
    if (EPV_observed < 10) {
      cat('❌ CRITICAL: EPV < 10 in Cox model\n')
      cat('   HR estimators are very unstable\n')
      cat('   Recommendation: Reduce number of variables or consider\n')
      cat('   penalization techniques (ridge Cox, lasso Cox)\n\n')
    } else if (EPV_observed < 20) {
      cat('⚠️  WARNING: EPV < 20 in Cox model\n')
      cat('   Cautious interpretation required\n')
      cat('   Recommendation: Report as limitation and perform\n')
      cat('   sensitivity analysis. Carefully verify proportional\n')
      cat('   hazards assumption\n\n')
    } else if (EPV_observed >= 20 && EPV_observed < 30) {
      cat('✓ Acceptable: EPV >= 20\n')
      cat('  Model meets minimum recommended for Cox\n\n')
    } else {
      cat('✓✓ Excellent: EPV >= 30\n')
      cat('   Cox model has robust estimates\n\n')
    }
    
    # Additional warning about censoring
    censoring_proportion <- n_censored / n_final
    if (censoring_proportion > 0.50) {
      cat('ℹ️  Note: High proportion of censoring (>', round(censoring_proportion*100), '%)\n')
      cat('   Verify if censoring is informative\n\n')
    }
  }
  
  return(invisible(list(
    n_final = n_final,
    events = n_events,
    censored = n_censored,
    k = k,
    EPV = EPV_observed,
    event_rate = observed_event_rate,
    language = language
  )))
}

#' @export
print.SampleSurvival <- function(x, ...) {
  if (x$language == 'es') {
    cat('\nTAMAÑO MUESTRAL TOTAL:', x$n_total, '\n')
    cat('  Eventos esperados:', x$events_needed, '\n')
    cat('  Censurados esperados:', x$censored_expected, '\n')
  } else {
    cat('\nTOTAL SAMPLE SIZE:', x$n_total, '\n')
    cat('  Events expected:', x$events_needed, '\n')
    cat('  Censored expected:', x$censored_expected, '\n')
  }
  invisible(x)
}

#' @export
print.SurvivalLogistics <- function(x, ...) {
  if (x$language == 'es') {
    cat('\nCohorte inicial:', x$initial_cohort, '\n')
    cat('Cohorte final:', x$final_cohort, '\n')
    cat('Duracion total:', x$total_duration_months, 'meses\n')
  } else {
    cat('\nInitial cohort:', x$initial_cohort, '\n')
    cat('Final cohort:', x$final_cohort, '\n')
    cat('Total duration:', x$total_duration_months, 'months\n')
  }
  invisible(x)
}

