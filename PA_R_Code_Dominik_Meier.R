library(quantmod)
library(tseries) 
library(PerformanceAnalytics) 
library(xts) 
library(Matrix)

tics <- c('^STOXX50E', 'GC=F', '^GSPC', 'CSBGC0.SW', 'IEF') 
e <- new.env()   # Erzeugen einer neuen Umgebung für die Daten
getSymbols(tics, from='2008-01-01', env = e)    # Daten herunterladen
TICS <- do.call(merge, eapply(env = e, Ad))    # Daten in xts-Objekt zusammenführen

# Samstage und Sonntage aus dem Datensatz entfernen
TICS <- TICS[!weekdays(index(TICS)) %in% c('Saturday', 'Sunday')]

# Fehlende Daten mit dem letzten gültigen Wert füllen
TICS <- na.locf(TICS)

# Löschen von Tagen mit komplett fehlenden Daten
TICS <- na.omit(TICS)

# Spaltennamen setzen
colnames(TICS) <- c('S&P 500', 'Swiss Bond', 'Gold', 'Eurostoxx50', 'US-Bond')
# Tägliche logarithmische Renditen berechnen
log_returns <- log(TICS / lag.xts(TICS, k = 1))

dates <- index(log_returns)
dates <- index(log_returns[-nrow(log_returns), ]) # -1 wegen NULL zeile
log_returns <- na.omit(log_returns)   # erste Zeile (NA) löschen
# Berechnung der einfachen täglichen Renditen
simple_returns <- (TICS / lag.xts(TICS, k = 1)) - 1
simple_returns <- na.omit(simple_returns)  # Erste Zeile (NA) löschen, da kein vorheriger Wert für den ersten Eintrag existiert

# Überprüfung 
head(simple_returns)
head(log_returns)
tail(log_returns)
###########################################################################

# Verschiedene Lambdas berechnet
# Halbwertszeiten für die Experte:  (Volas und Korrelationen)
# 10/21; 21/63; 63/125; 125/250; 250/500

compute_lambda <- function(halflife) {
  exp(log(0.5) / halflife)
}

halflives <- c(10, 21, 63, 125, 250, 500)

for (hl in halflives) {
  assign(paste0("lambda", hl), compute_lambda(hl))
}
rm(hl)

#########################################################################

# EWMA Kovarianz schätzen

ewma_covariance_all <- function(log_returns, lambda_vol, lambda_corr) {
  n <- ncol(log_returns)
  T <- nrow(log_returns)
  cov_matrices_list <- list()
  vol_estimates <- matrix(1, n, T)  # Volatilitäten initialisieren
  S <- cov(log_returns)  # Anfängliche Kovarianzschätzung
  
  for (t in 1:(T-1)) {
    rt <- as.matrix(log_returns[t, ])
    
    # Schritt 1: Volatilitäten mit EWMA schätzen
    vol_estimates[, t + 1] <- sqrt(lambda_vol * (vol_estimates[, t]^2) + (1 - lambda_vol) * rt^2)
    
    # Schritt 2: Renditen standardisieren
    rt_standardized <- sweep(rt, 2, vol_estimates[, t], FUN = "/")
    
    # Schritt 3: EWMA-Schätzung der Kovarianz der standardisierten Renditen
    for (i in 1:n) {
      for (j in i:n) {
        S[i, j] <- lambda_corr * S[i, j] + (1 - lambda_corr) * rt_standardized[i] * rt_standardized[j]
        if (i != j) {
          S[j, i] <- S[i, j]
        }
      }
    }
    
    # Schritt 4: Skalierung zur Bildung der Korrelationsmatrix
    S_scaled <- S / sqrt(outer(diag(S), diag(S), FUN = "*"))
    
    # Schritt 5: Endgültigen Kovarianzschätzer bilden
    final_cov <- vol_estimates[, t] %*% t(vol_estimates[, t]) * S_scaled
    
    cov_matrices_list[[length(cov_matrices_list) + 1]] <- final_cov
  }
  
  return(cov_matrices_list)
}



# Halbwertszeiten für die Lambda-Werte
halflife_vol <- c(10, 21, 63, 125, 250)
halflife_corr <- c(21, 63, 125, 250, 500)

lambda_vol_values <- c(lambda10, lambda21, lambda63, lambda125, lambda250)
lambda_corr_values <- c(lambda21, lambda63, lambda125, lambda250, lambda500)

# Liste der gewünschten Paare von Halbwertszeiten
desired_pairs <- list(c(10,21), c(21,63), c(63,125), c(125,250), c(250,500))

# Liste, um die Ergebnisse zu speichern
results_all <- list()

for (pair in desired_pairs) {
  i <- which(halflife_vol == pair[1])
  j <- which(halflife_corr == pair[2])
  lambda_vol_value <- lambda_vol_values[i]
  lambda_corr_value <- lambda_corr_values[j]
  key <- paste0("lambda_vol", pair[1], "_lambda_corr", pair[2])
  
  cov_matrices_list <- ewma_covariance_all(log_returns, lambda_vol_value, lambda_corr_value)
  results_all[[key]] <- cov_matrices_list
}

# Funktion zum Überprüfen, ob eine Matrix ausschließlich Nullen enthält
is_all_zero <- function(matrix) {
  all(matrix == 0)
}

# Kovarianzmatrizen anpassen
adjust_cov_matrices <- function(cov_list) {
  len <- length(cov_list)
  if (len > 1 && is_all_zero(cov_list[[len]])) {
    cov_list[[len]] <- cov_list[[len - 1]]  # Ersetze die letzte Matrix durch die vorletzte, falls nötig
  }
  
  return(cov_list)
}

# Anwenden der Funktion auf jede Liste in 'results_all'
results_all <- lapply(results_all, adjust_cov_matrices)

# Zuweisung neuer Namen zu den Listen, falls erforderlich
names(results_all) <- c("cov_ewma_1", "cov_ewma_2", "cov_ewma_3", "cov_ewma_4", "cov_ewma_5")


# Übertragen der aktualisierten Listen in die globale Umgebung
list2env(results_all, envir = .GlobalEnv)

head(cov_ewma_1)
tail(cov_ewma_5)

###########################################################################
########*************** Positiv semi definit machen eigenwert Anpassung *******************

# Laden des Matrix-Pakets für die Funktion nearPD zur Bereinigung der Kovarianzmatrix
library(Matrix)

# Funktion, um sicherzustellen, dass die Kovarianzmatrix symmetrisch ist
make_symmetric <- function(cov_matrix) {
  symmetric_cov_matrix <- 0.5 * (cov_matrix + t(cov_matrix))
  return(symmetric_cov_matrix)
}

# Funktion, um die Kovarianzmatrix positiv definit zu machen
make_positive_definite <- function(cov_matrix) {
  cleaned_cov_matrix <- nearPD(cov_matrix)$mat
  return(cleaned_cov_matrix)
}

# Anwenden der Funktion auf jede Kovarianzmatrix in jeder Liste
results_all <- lapply(results_all, function(cov_list) {
  lapply(cov_list, function(cov_matrix) {
    # Sicherstellen, dass die Kovarianzmatrix symmetrisch ist
    symmetric_cov_matrix <- make_symmetric(cov_matrix)
    
    # Kovarianzmatrix in positive definite Form umwandeln
    positive_definite_cov_matrix <- make_positive_definite(symmetric_cov_matrix)
    
    return(positive_definite_cov_matrix)
  })
})

# Übertragen der aktualisierten Listen in die globale Umgebung
list2env(results_all, envir = .GlobalEnv)

# Beispiel, um die erste Kovarianzmatrix von cov_ewma_1 zu überprüfen
head(cov_ewma_1[[1]])




###************Funktion zum Überprüfen der positiven Semidefinitheit**********************************

is_positive_semidefinite <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  all(eigenvalues >= 0L)
}

# Überprüfen jeder Kovarianzmatrix in der Liste cov_ewma_1 auf positive Semidefinitheit
psd_checks_ewma_1 <- sapply(cov_ewma_1, is_positive_semidefinite)

# Ausgabe der Ergebnisse
print(psd_checks_ewma_1)

# Funktion, die die Summe der negativen Eigenwerte einer Matrix berechnet
sum_negative_eigenvalues <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  sum_negatives <- sum(eigenvalues[eigenvalues < 0])
  return(sum_negatives)
}

# Anwenden der Funktion auf jede Kovarianzmatrix in cov_ewma_1
neg_eigenvalues_sums <- sapply(cov_ewma_1, sum_negative_eigenvalues)

# Ausgabe der Summe der negativen Eigenwerte für jede Kovarianzmatrix
print(neg_eigenvalues_sums)
min(neg_eigenvalues_sums)

##********************  Gewichte für MVP berechnen *****************************

# Laden des NMOF-Pakets für die MVP-Berechnung
library(NMOF)

# Funktion zur Anpassung der Gewichte
adjust_weights <- function(weights) {
  weights[weights < 0] <- 0  # Setze alle negativen Gewichte auf Null
  weights / sum(weights)     # Normiere die Gewichte, so dass ihre Summe 1 ergibt
}

# Listen, um die angepassten MVP-Gewichte und Varianzen separat zu speichern
mvp_weights_list <- list()
mvp_variances_list <- list()

# Iteration durch jede Liste von täglichen Kovarianzmatrizen
for(i in 1:5) {
  cov_list_name <- paste0("cov_ewma_", i)
  cov_list <- get(cov_list_name)  
  
  weights_list <- list()     # Liste für Gewichte
  variances_list <- list()   # Liste für Varianzen
  
  # Iteration durch jede Kovarianzmatrix in der aktuellen Liste
  for(j in seq_along(cov_list)) {
    # Berechnung der MVP-Gewichte und Varianz für die aktuelle Kovarianzmatrix
    mvp_result <- minvar(cov_list[[j]], wmin = 0, wmax = 1)###########################################
    #S&P 500', 'Swiss Bond', 'Gold', 'Eurostoxx50', 'US-Bond'
    # Anpassung der Gewichte, um negative Gewichte zu vermeiden
    adjusted_weights <- adjust_weights(mvp_result)
    
    # Speichern der angepassten Gewichte und Varianz in separaten Listen
    weights_list[[j]] <- adjusted_weights
    variances_list[[j]] <- attr(mvp_result, "variance")
  }
  
  # Speichern der Listen der Gewichte und Varianzen in den Hauptlisten
  mvp_weights_list[[cov_list_name]] <- weights_list
  mvp_variances_list[[cov_list_name]] <- variances_list
}

#Listen umbennenen
for(i in 1:5) {
  old_weights_list_name <- paste0("cov_ewma_", i)
  new_weights_list_name <- paste0("mvp_w_", i)
  assign(new_weights_list_name, mvp_weights_list[[old_weights_list_name]])
  
  old_variances_list_name <- paste0("cov_ewma_", i)
  new_variances_list_name <- paste0("mvp_v_", i)
  assign(new_variances_list_name, mvp_variances_list[[old_weights_list_name]])
}

##************** Gewichte als XTS objekt abspeichern ************************


# Definieren der Spaltennamen
asset_names <- c('S&P 500', 'Swiss Bond', 'Gold', 'Eurostoxx50', 'US-Bond')

# Erstellen einer Funktion zur Konvertierung einer Liste von Gewichten in ein XTS-Objekt
# und Benennung der Spalten
convert_to_xts_and_name <- function(weights_list, dates, col_names) {
  # Umwandeln der Liste in eine Matrix
  weights_matrix <- do.call(rbind, weights_list)
  
  # Erstellen des XTS-Objekts
  weights_xts <- xts(weights_matrix, order.by = dates)
  
  # Benennen der Spalten
  colnames(weights_xts) <- col_names
  
  return(weights_xts)
}
dates <- index(log_returns[-nrow(log_returns), ]) # -1 wegen NULL zeile
# Anwendung der Funktion auf jede MVP-Gewichtsliste
xts_mvp_w_1 <- convert_to_xts_and_name(mvp_w_1, dates, asset_names)
xts_mvp_w_2 <- convert_to_xts_and_name(mvp_w_2, dates, asset_names)
xts_mvp_w_3 <- convert_to_xts_and_name(mvp_w_3, dates, asset_names)
xts_mvp_w_4 <- convert_to_xts_and_name(mvp_w_4, dates, asset_names)
xts_mvp_w_5 <- convert_to_xts_and_name(mvp_w_5, dates, asset_names)

# Überprüfung der ersten paar Zeilen eines der XTS-Objekte
tail(xts_mvp_w_5)
tail(xts_mvp_w_1)
##************** monatlichen Gewichte MVPs ****************************************

library(xts)

# Definieren der Spaltennamen
asset_names <- c('S&P 500', 'Swiss Bond', 'Gold', 'Eurostoxx50', 'US-Bond')

# Funktion zur Konvertierung einer Liste von Gewichten in ein XTS-Objekt und Benennung der Spalten
convert_to_xts_and_name <- function(weights_list, dates, col_names) {
  # Umwandeln der Liste in eine Matrix
  weights_matrix <- do.call(rbind, weights_list)
  
  # Erstellen des XTS-Objekts
  weights_xts <- xts(weights_matrix, order.by = dates)
  
  # Benennen der Spalten
  colnames(weights_xts) <- col_names
  
  # Extrahieren der Gewichte am Ende jeden Monats
  month_end_weights_xts <- weights_xts[endpoints(weights_xts, on = "months"), ]
  
  return(month_end_weights_xts)
}

# Anwendung der Funktion auf jede MVP-Gewichtsliste
m_w_1 <- convert_to_xts_and_name(mvp_w_1, dates, asset_names)
m_w_2 <- convert_to_xts_and_name(mvp_w_2, dates, asset_names)
m_w_3 <- convert_to_xts_and_name(mvp_w_3, dates, asset_names)
m_w_4 <- convert_to_xts_and_name(mvp_w_4, dates, asset_names)
m_w_5 <- convert_to_xts_and_name(mvp_w_5, dates, asset_names)

# Überprüfung der ersten paar Zeilen eines der XTS-Objekte
head(m_w_1)
tail(m_w_5)

par(mfrow = c(3, 2))
chart.StackedBar(
  m_w_2,
  colorset = NULL,
  space = 0.2,
  cex.axis = 0.8,
  cex.legend = 0.8,
  cex.lab = 1,
  cex.labels = 0.8,
  cex.main = 1,
  xaxis = TRUE,
  legend.loc = "under",
  element.color = "darkgray",
  unstacked = TRUE,
  xlab = "Date",
  ylab = "Value",
  ylim = NULL,
  date.format = "%b %y",
  major.ticks = "auto",
  minor.ticks = TRUE,
  las = 0,
  xaxis.labels = NULL,
)

##************** monatlichen kummulierte Renditen log_returns ************************

# Funktion zur Berechnung der kumulierten monatlichen Rendite aus Log-Renditen
cumulative_monthly_return <- function(column) {
  apply.monthly(column, sum)
}

# Berechnen der monatlichen Renditen für jedes Portfolio
monthly_returns_list <- lapply(log_returns, cumulative_monthly_return)

# Ergebnisse in ein xts-Objekt umwandeln
monthly_log_returns_xts <- do.call(merge, monthly_returns_list)

# Überprüfen der berechneten monatlichen Renditen
head(monthly_log_returns_xts)
tail(monthly_log_returns_xts)
tail(log_returns)

##************** MVP Renditen berechnen In-Sample daily ************************

# Funktion zur Berechnung der täglichen Renditen für ein MVP
calculate_daily_mvp_returns <- function(weights_xts, log_returns) {
  # Multiplizieren der log-returns mit den Gewichten für jedes Asset
  weighted_log_returns <- sweep(log_returns, 1, weights_xts, "*")
  
  # Aggregieren der gewichteten log-returns über alle Assets
  daily_mvp_returns <- rowSums(weighted_log_returns, na.rm = TRUE)
  
  return(daily_mvp_returns)
}

# Berechnung der täglichen Renditen für jedes MVP
daily_returns_mvp_1 <- calculate_daily_mvp_returns(xts_mvp_w_1, log_returns)
daily_returns_mvp_2 <- calculate_daily_mvp_returns(xts_mvp_w_2, log_returns)
daily_returns_mvp_3 <- calculate_daily_mvp_returns(xts_mvp_w_3, log_returns)
daily_returns_mvp_4 <- calculate_daily_mvp_returns(xts_mvp_w_4, log_returns)
daily_returns_mvp_5 <- calculate_daily_mvp_returns(xts_mvp_w_5, log_returns)

# Umwandlung der täglichen Renditen in xts-Objekte
xts_daily_returns_mvp_1 <- xts(daily_returns_mvp_1, order.by = index(log_returns))
xts_daily_returns_mvp_2 <- xts(daily_returns_mvp_2, order.by = index(log_returns))
xts_daily_returns_mvp_3 <- xts(daily_returns_mvp_3, order.by = index(log_returns))
xts_daily_returns_mvp_4 <- xts(daily_returns_mvp_4, order.by = index(log_returns))
xts_daily_returns_mvp_5 <- xts(daily_returns_mvp_5, order.by = index(log_returns))

# Überprüfung der ersten paar Zeilen eines der XTS-Objekte
head(xts_daily_returns_mvp_1)


merged_mvp_xts <- merge.xts(xts_daily_returns_mvp_1, xts_daily_returns_mvp_2,
                            xts_daily_returns_mvp_3 , xts_daily_returns_mvp_4,
                            xts_daily_returns_mvp_5)

# Umbenennen der Spalten
colnames(merged_mvp_xts) <- c("10/21", "21/63", "63/125", "125/250", "250/500")

# Überprüfen des kombinierten xts-Objekts
head(merged_mvp_xts)
tail(merged_mvp_xts)


# Visualisierung IN-Sample daily rebalancing
library(RColorBrewer)
colorset1 <- brewer.pal(n = 5, name = "Set1")
charts.PerformanceSummary(merged_mvp_xts, colorset = colorset1)


##************** IN-Sample Monthly mit log_returns ************************


# Funktion zur Berechnung der monatlichen Renditen für ein MVP
calculate_monthly_mvp_returns <- function(weights_xts, returns_xts) {
  mvp_monthly_returns <- NULL
  
  len <- min(nrow(weights_xts), nrow(returns_xts))
  
  # Iterieren durch die Monate
  for (i in 1:len) {
    # Multiplizieren der Gewichte mit den Renditen und Summieren
    monthly_return <- sum(returns_xts[i, ] * weights_xts[i, ])
    mvp_monthly_returns <- c(mvp_monthly_returns, monthly_return)
  }
  
  return(xts(mvp_monthly_returns, order.by = index(returns_xts)[1:length(mvp_monthly_returns)]))
}

# Berechnen der monatlichen Renditen für jedes MVP
monthly_returns_mvp_list <- list()
for (i in 1:5) {
  weights_xts <- get(paste0("m_w_", i))
  monthly_returns_mvp_list[[i]] <- calculate_monthly_mvp_returns(weights_xts, monthly_log_returns_xts)
}

# Zusammenführen der monatlichen Renditen aller MVPs in einem XTS-Objekt
all_mvp_returns <- do.call(merge, monthly_returns_mvp_list)
colnames(all_mvp_returns) <- paste0("MVP_", 1:5)

# Überprüfung der monatlichen Renditen
head(all_mvp_returns)
tail(all_mvp_returns)

# Farbset für die Diagramme festlegen
colorset1 <- brewer.pal(n = 5, name = "Set1")

# Erstellen der Performance-Charts
charts.PerformanceSummary(all_mvp_returns, colorset = colorset1)


start_date1 <- as.Date("2010-01-01")
IN_Portfolio <- window(all_mvp_returns, start=start_date1)
colnames(IN_Portfolio) <- c("10/21", "21/63", "63/125", "125/250", "250/500")
charts.PerformanceSummary(IN_Portfolio, colorset = colorset1)

##******* Visualisierung der monatlichen Renditen und Anpassung der Liste ************

# Berechnung der kumulativen Renditen
cumulative_returns <- apply(simple_returns + 1, 2, cumprod) - 1

library(reshape2)
library(ggplot2)
# Umwandeln der Daten in ein langes Format für ggplot



############****** MVP tägliche Renditen XTS Zeitreihen berechenen *********************************
log_returns <- log_returns[-nrow(log_returns), ] # -1 wegen NULL zeile
# Berechnung der täglichen Renditen für jedes MVP
calculate_daily_mvp_returns <- function(log_returns, mvp_weights) {
  daily_mvp_returns <- rep(NA, nrow(log_returns))
  
  for (i in 1:length(mvp_weights)) {
    # Multiplikation der täglichen Renditen mit den entsprechenden Gewichten
    daily_mvp_returns[i] <- sum(log_returns[i, ] * mvp_weights[[i]])
  }
  
  return(daily_mvp_returns)
}

# Umwandeln der MVP-Renditen in xts-Objekte
create_mvp_xts <- function(daily_returns, dates) {
  xts_obj <- xts(daily_returns, order.by = as.Date(dates))
  return(xts_obj)
}

# Berechnung der täglichen Renditen für jedes MVP und Speicherung als xts-Objekt
mvp_xts_list <- list()

for(i in 1:5) {
  mvp_weights <- get(paste0("mvp_w_", i))
  daily_mvp_returns <- calculate_daily_mvp_returns(log_returns, mvp_weights)
  
  # Umwandlung in xts
  mvp_xts <- create_mvp_xts(daily_mvp_returns, dates)
  mvp_xts_list[[paste0("mvp_xts_", i)]] <- mvp_xts
}

# Zuweisung der xts-Objekte in die globale Umgebung
list2env(mvp_xts_list, envir = .GlobalEnv)

# Zusammenführen der einzelnen MVP xts-Objekte
merged_mvp_xts <- merge.xts(mvp_xts_1, mvp_xts_2, mvp_xts_3, mvp_xts_4, mvp_xts_5)

# Umbenennen der Spalten
colnames(merged_mvp_xts) <- c("10/21", "21/63", "63/125", "125/250", "250/500")

# Überprüfen des kombinierten xts-Objekts
head(merged_mvp_xts)
tail(merged_mvp_xts)

############****** Rollbacktesting jeden Monat Rebalancing ab 2008 ******************************************

# Funktion zur Berechnung der monatlichen Renditen für ein MVP mit Verzögerung
calculate_monthly_mvp_returns_lagged <- function(weights_xts, returns_xts) {
  # Gewichte um einen Monat verzögern
  lagged_weights_xts <- lag.xts(weights_xts, k = +1)
  
  # Initialisieren der Matrix für monatliche Renditen
  mvp_monthly_returns <- matrix(nrow = nrow(lagged_weights_xts), ncol = 1)
  
  for (i in 2:nrow(lagged_weights_xts)) {  
    # Monatliche Renditen berechnen
    mvp_monthly_returns[i, 1] <- sum(returns_xts[i, ] * lagged_weights_xts[i, ])
  }
  
  return(xts(mvp_monthly_returns, order.by = index(returns_xts)))
}

# Anwendung auf alle MVPs
all_mvp_monthly_returns_lagged <- list()
for (i in 1:5) {
  weights_xts <- get(paste0("m_w_", i))
  all_mvp_monthly_returns_lagged[[i]] <- calculate_monthly_mvp_returns_lagged(weights_xts, monthly_log_returns_xts)
}

# Zusammenführen der monatlichen Renditen aller MVPs in einem XTS-Objekt
all_mvp_returns_lagged_xts <- do.call(merge, all_mvp_monthly_returns_lagged)
colnames(all_mvp_returns_lagged_xts) <- paste0("MVP_", 1:5)

# Überprüfung der monatlichen Renditen
head(all_mvp_returns_lagged_xts)
tail(all_mvp_returns_lagged_xts)


charts.PerformanceSummary(all_mvp_returns_lagged_xts, colorset = colorset1)


############****** Rollbacktesting jeden Monat Rebalancing ab 2010 ***********************************

start_date <- as.Date("2010-01-01")
all_mvp_returns_lagged_xts_2010_onwards <- window(all_mvp_returns_lagged_xts, start=start_date)

# Überprüfung der gefilterten monatlichen Renditen
head(all_mvp_returns_lagged_xts_2010_onwards)
tail(all_mvp_returns_lagged_xts_2010_onwards)

# Visualisierung der Performance ab 2010
colorset1 <- brewer.pal(n = 5, name = "Set1")
charts.PerformanceSummary(all_mvp_returns_lagged_xts_2010_onwards, colorset = colorset1)

############****** Analyse MVPs *************************************

MVP_Comp <- all_mvp_returns_lagged_xts_2010_onwards

start_date1 <- as.Date("2010-01-01")
MVP_Comp <- window(MVP_Comp , start=start_date1)


chart.Drawdown(
  MVP_Comp,
  geometric = TRUE,
  legend.loc = "bottomleft",
  colorset = (2:7),
  plot.engine = "default",
  main= "Max Drawdown MVPs" # nicht benutzen
)



charts.TimeSeries(allop, legend.loc = "bottomleft")

VaR(
  R = MVP_Comp,
  p = 0.95,
  method = c("modified", "gaussian", "historical", "kernel"),
  clean = c("none", "boudt", "geltner", "locScaleRob"),
  portfolio_method = c("single"),
  weights = NULL,
  mu = NULL,
  sigma = NULL,
  m3 = NULL,
  m4 = NULL,
  invert = TRUE,
  SE = FALSE,
  SE.control = NULL)

DownsideDeviation(
  MVP_Comp,
  MAR = 0,
  method = c("full", "subset"),
  potential = FALSE
)


VolatilitySkewness(MVP_Comp, MAR = 0, stat = c("volatility", "variability"))

m_w_1 <- window(m_w_1, start=start_date1)
library(RColorBrewer)
colnames(monthly_mvp_weights_xts) <- c('S&P 500', 'Swiss Bond', 'Gold', 'Eurostoxx50', 'US-Bond')
chart.StackedBar(
  monthly_mvp_weights_xts,
  main= "Portfolio Weights MVP6 (Sample Covariance)",
  colorset = c("red",4,"lightgreen","darkblue","orange"),
  space = 0.2,
  cex.axis = 0.8,
  cex.legend = 0.8,
  cex.lab = 1,
  cex.labels = 0.8,
  cex.main = 1,
  xaxis = TRUE,
  legend.loc = "under",
  element.color = "darkgray",
  unstacked = TRUE,
  xlab = "Date",
  ylab = "Value",
  ylim = NULL,
  date.format = "%Y",
  major.ticks = "auto",
  minor.ticks = TRUE,
  las = 0,
  xaxis.labels = NULL,
)
OmegaSharpeRatio(allop)

chart.TimeSeries(allop)
chart.StackedBar(
  m_w_5,
  colorset = NULL,
  space = 0.2,
  cex.axis = 0.8,
  cex.legend = 0.8,
  cex.lab = 1,
  cex.labels = 0.8,
  cex.main = 1,
  xaxis = TRUE,
  legend.loc = "under",
  element.color = "darkgray",
  unstacked = TRUE,
  xlab = "Date",
  ylab = "Value",
  ylim = NULL,
  date.format = "%b %y",
  major.ticks = "auto",
  minor.ticks = TRUE,
  las = 0,
  xaxis.labels = NULL,
)


chart.RiskReturnScatter(
  MVP_Comp,
  main = "Annualized Return and Risk",
  xlab = "Annualized Risk",
  ylab = "Annualized Return",
  method = "calc",
  geometric = TRUE,
  scale = NA,
)

  
chart.VaRSensitivity(
  MVP_Comp)
  

StdDev(MVP_Comp)

a_barplot <- StdDev.annualized(MVP_Comp, scale = NA)
barplot(a_barplot)

MSquared() # Hier vergleichen min risk der MVP gegen sample covariance
# M squared excess is the quantity above the standard M. There is a geometric excess return which is
# better for Bacon and an arithmetic excess return

DrawdownDeviation(MVP_Comp)
drawdownplot <- DrawdownDeviation(MVP_Comp)
barplot(drawdownplot)

#annualisiert
chart.Bar(MVP_Comp, legend.loc = "bottom", colorset = (1:5))

chart.Boxplot(MVP_Comp)



# Calculates Expected Shortfall(ES)
ETL_op <- ETL(
  R = OP,
  p = 0.95,
  method = c("modified", "gaussian", "historical"),
  clean = c("none", "boudt", "geltner", "locScaleRob"),
  portfolio_method = c("single", "component"),
  weights = NULL,
  mu = NULL,
  sigma = NULL,
  m3 = NULL,
  m4 = NULL,
  invert = TRUE,
  operational = TRUE,
  SE = FALSE,
  SE.control = NULL
)

round(ETL_op,4)
############****** sample covariance monatlich **********************************
# Funktion zum Berechnen der monatlichen MVP-Gewichte

#log_returns <- log_returns[-nrow(log_returns), ] # falls 1. von monat
#log_returns <- log_returns[-nrow(log_returns-1), ] # falls 2. von monat
#log_returns <- log_returns[-nrow(log_returns-2), ] # falls 3. von monat
calculate_monthly_mvp_weights <- function(log_returns) {
  endpoints <- endpoints(log_returns, on = "months")
  monthly_weights <- list()
  
  for (i in 1:(length(endpoints) - 1)) {
    # Extrahieren der Log-Renditen für den aktuellen Monat
    monthly_returns <- log_returns[(endpoints[i] + 1):endpoints[i + 1], ]
    
    # Berechnen der Kovarianzmatrix für den aktuellen Monat
    sample_cov_matrix <- cov(monthly_returns)
    
    # Berechnen der MVP-Gewichte
    mvp_weights <- minvar(sample_cov_matrix, wmin = 0, wmax = 1)
    
    monthly_weights[[i]] <- mvp_weights
  }
  
  return(do.call(rbind, monthly_weights))
}

# Berechnen der monatlichen MVP-Gewichte
monthly_mvp_weights <- calculate_monthly_mvp_weights(log_returns)

# Erstellen eines xts-Objekts mit den monatlichen Gewichten
monthly_mvp_weights_xts <- xts(monthly_mvp_weights, order.by = index(log_returns)[endpoints(log_returns, on = "months")[-1]])
colnames(monthly_mvp_weights_xts) <- asset_names


# Funktion zur Berechnung der monatlichen Renditen des MVP mit verzögerten Gewichten
calculate_monthly_mvp_returns_lagged <- function(weights_xts, log_returns_xts) {
  # Gewichte um einen Monat verzögern
  lagged_weights_xts <- lag.xts(monthly_mvp_weights_xts, k = +1)
  
  monthly_mvp_returns <- NULL
  
  # Gruppierung der Log-Renditen nach Monaten
  monthly_log_returns <- apply.monthly(log_returns, colSums)
  
  for (i in 2:nrow(lagged_weights_xts)) {
    # Extrahieren der verzögerten Gewichte
    lagged_weights <- as.numeric(lagged_weights_xts[i, ])
    
    # Extrahieren der monatlichen summierten Log-Renditen für den aktuellen Monat
    if (index(monthly_log_returns)[i] %in% index(log_returns)) {
      current_monthly_log_returns <- monthly_log_returns[index(lagged_weights_xts)[i], ]
      
      # Berechnen der monatlichen Rendite für das Portfolio
      monthly_return <- sum(current_monthly_log_returns * lagged_weights, na.rm = TRUE)
      monthly_mvp_returns <- c(monthly_mvp_returns, monthly_return)
    }
  }
  
  return(xts(monthly_mvp_returns, order.by = index(lagged_weights_xts)[-1]))
}


# Berechnen der monatlichen Renditen des MVP mit verzögerten Gewichten
monthly_mvp_returns_lagged <- calculate_monthly_mvp_returns_lagged(monthly_mvp_weights_xts, log_returns)
mvp_w_m <- monthly_mvp_weights_xts
# Aktualisieren des MVP-Vergleichsobjekts
MVP_Comp_updated <- merge.xts(MVP_Comp, monthly_mvp_returns_lagged)
colnames(MVP_Comp_updated)[ncol(MVP_Comp_updated)] <- "Sample_Cov_m"

# Visualisieren der Portfolio-Leistung
charts.PerformanceSummary(MVP_Comp_updated, colorset = 2:7)



############****** sample covariance historisch **********************************

calculate_monthly_mvp_weights <- function(log_returns) {
  endpoints <- endpoints(log_returns, on = "months")
  monthly_weights <- list()
  
  for (i in 1:(length(endpoints) - 1)) {
    # Extrahieren der Log-Renditen bis zum aktuellen Monat
    cumulative_returns <- log_returns[1:endpoints[i + 1], ]
    
    # Berechnen der kumulativen Kovarianzmatrix bis zum aktuellen Monat
    cumulative_cov_matrix <- cov(cumulative_returns)
    
    # Berechnen der MVP-Gewichte
    mvp_weights <- minvar(cumulative_cov_matrix, wmin = 0, wmax = 1)
    
    monthly_weights[[i]] <- mvp_weights
  }
  
  return(do.call(rbind, monthly_weights))
}


# Berechnen der monatlichen MVP-Gewichte
monthly_mvp_weights <- calculate_monthly_mvp_weights(log_returns)

# Erstellen eines xts-Objekts mit den monatlichen Gewichten
monthly_mvp_weights_xts <- xts(monthly_mvp_weights, order.by = index(log_returns)[endpoints(log_returns, on = "months")[-1]])
colnames(monthly_mvp_weights_xts) <- asset_names


# Funktion zur Berechnung der monatlichen Renditen des MVP mit verzögerten Gewichten
calculate_monthly_mvp_returns_lagged <- function(weights_xts, log_returns_xts) {
  # Gewichte um einen Monat verzögern
  lagged_weights_xts <- lag.xts(monthly_mvp_weights_xts, k = +1)
  
  monthly_mvp_returns <- NULL
  
  # Gruppierung der Log-Renditen nach Monaten
  monthly_log_returns <- apply.monthly(log_returns, colSums)
  
  for (i in 2:nrow(lagged_weights_xts)) {
    # Extrahieren der verzögerten Gewichte
    lagged_weights <- as.numeric(lagged_weights_xts[i, ])
    
    # Extrahieren der monatlichen summierten Log-Renditen für den aktuellen Monat
    if (index(monthly_log_returns)[i] %in% index(log_returns)) {
      current_monthly_log_returns <- monthly_log_returns[index(lagged_weights_xts)[i], ]
      
      # Berechnen der monatlichen Rendite für das Portfolio
      monthly_return <- sum(current_monthly_log_returns * lagged_weights, na.rm = TRUE)
      monthly_mvp_returns <- c(monthly_mvp_returns, monthly_return)
    }
  }
  
  return(xts(monthly_mvp_returns, order.by = index(lagged_weights_xts)[-1]))
}


# Berechnen der monatlichen Renditen des MVP mit verzögerten Gewichten
monthly_mvp_returns_lagged_h <- calculate_monthly_mvp_returns_lagged(monthly_mvp_weights_xts, log_returns)
mvp_w_h <- monthly_mvp_weights_xts
# Aktualisieren des MVP-Vergleichsobjekts
MVP_Comp_updated <- merge.xts(MVP_Comp_updated, monthly_mvp_returns_lagged_h)
colnames(MVP_Comp_updated)[ncol(MVP_Comp_updated)] <- "Sample_Cov_h"

# Visualisieren der Portfolio-Leistung
charts.PerformanceSummary(MVP_Comp_updated, colorset = 2:8)

colnames(MVP_Comp) <- c("10/21", "21/63", "63/125", "125/250", "250/500")
colnames(MVP_Comp_updated) <- c("10/21", "21/63", "63/125", "125/250", "250/500", "Sample_Cov")
 ############****** Analyse MVPs Sample *************************************

# lösche alle Zeilen die NA sind
na.omit(MVP_Comp_updated)

start_date1 <- as.Date("2010-01-01")
MVP_Comp_updated1 <- window(MVP_Comp_updated1 , start=start_date1)

charts.RollingPerformance(
  MVP_Comp_updated,
  width = 12,
  Rf = 0,
  main = "Rolling 12 month Performance Chart von allen MVPs",
  event.labels = FALSE,
  legend.loc = "bottomleft")

chart.Drawdown(
  IN_Portfolio,
  geometric = TRUE,
  legend.loc = FALSE,
  colorset = colors,
  plot.engine = "default",
  main= "Max Drawdown MVPs" # nicht benutzen
)

chart.CumReturns(
  IN_Portfolio,
  main = "Kummulierte Log-Renditen der MVPs (in-sample)",
  wealth.index = FALSE,
  geometric = TRUE,
  legend.loc = "topleft",
  ylab = "Log-Rendite",
  colorset = colors,
  begin = c("first", "axis"),
  plot.engine = "default"
)



VaR(
  R = OP,
  p = 0.95,
  method = c("modified", "gaussian", "historical", "kernel"),
  clean = c("none", "boudt", "geltner", "locScaleRob"),
  portfolio_method = c("single"),
  weights = NULL,
  mu = NULL,
  sigma = NULL,
  m3 = NULL,
  m4 = NULL,
  invert = TRUE,
  SE = FALSE,
  SE.control = NULL)

DownsideDeviation(
  MVP_Comp_updated,
  MAR = 0,
  method = c("full", "subset"),
  potential = FALSE
)


VolatilitySkewness(MVP_Comp_updated, MAR = 0, stat = c("volatility", "variability"))




chart.StackedBar(
  m_w_4,
  colorset = 10:15,
  space = 0.2,
  cex.axis = 0.8,
  cex.legend = 0.8,
  cex.lab = 1,
  cex.labels = 0.8,
  cex.main = 1,
  xaxis = TRUE,
  legend.loc = "under",
  element.color = "darkgray",
  unstacked = TRUE,
  xlab = "Date",
  ylab = "Value",
  ylim = NULL,
  date.format = "%b %y",
  major.ticks = "auto",
  minor.ticks = TRUE,
  las = 0,
  xaxis.labels = NULL,
)

chart.StackedBar(
  mvp_w_dcc,
  colorset = NULL,
  space = 0.2,
  cex.axis = 0.8,
  cex.legend = 0.8,
  cex.lab = 1,
  cex.labels = 0.8,
  cex.main = 1,
  xaxis = TRUE,
  legend.loc = "under",
  element.color = "darkgray",
  unstacked = TRUE,
  xlab = "Date",
  ylab = "Value",
  ylim = NULL,
  date.format = "%b %y",
  major.ticks = "auto",
  minor.ticks = TRUE,
  las = 0,
  xaxis.labels = NULL,
)


chart.RiskReturnScatter(
  MVP_Comp_updated,
  main = "Annualized Return and Risk",
  xlab = "Annualized Risk",
  ylab = "Annualized Return",
  method = "calc",
  geometric = TRUE,
  scale = NA,
)


chart.VaRSensitivity(
  MVP_Comp_updated)


StdDev(MVP_Comp_updated1)

a_barplot <- StdDev.annualized(MVP_Comp_updated, scale = NA)
barplot(a_barplot)

MSquared() # Hier vergleichen min risk der MVP gegen sample covariance
# M squared excess is the quantity above the standard M. There is a geometric excess return which is
# better for Bacon and an arithmetic excess return

DrawdownDeviation(MVP_Comp_updated)
drawdownplot <- DrawdownDeviation(MVP_Comp_updated)
barplot(drawdownplot)

#annualisiert
chart.Bar(MVP_Comp_updated1, colorset = colors, main = "Monatliche Volatilität MVPs ")

chart.Boxplot(MVP_Comp_updated)



# Calculates Expected Shortfall(ES)
ETL(
  R = MVP_Comp_updated,
  p = 0.95,
  method = c("modified", "gaussian", "historical"),
  clean = c("none", "boudt", "geltner", "locScaleRob"),
  portfolio_method = c("single", "component"),
  weights = NULL,
  mu = NULL,
  sigma = NULL,
  m3 = NULL,
  m4 = NULL,
  invert = TRUE,
  operational = TRUE,
  SE = FALSE,
  SE.control = NULL
)

############****** Analyse IST und SOLL*************************************
library(zoo)
library(highcharter)


all_mvp_returns_lagged_xts <- na.locf(all_mvp_returns_lagged_xts)

plot(MVP_Comp_updated)

window <- 12

pr_rolling_sd_24 <- rollapply(MVP_Comp_updated$`Sample Cov`, FUN = sd,
                              width = window)  

head(pr_rolling_sd_24)

highchart(type = "stock") %>% 
  hc_title(text = "12-Month Rolling volatility") %>% 
  hc_add_series(pr_rolling_sd_24) %>% 
  hc_add_theme(hc_theme_flat()) %>% 
  hc_yAxis(
    labels = list(format = "{value}%"), 
    opposite = FALSE) %>%
  hc_navigator(enabled = FALSE) %>% 
  hc_scrollbar(enabled = FALSE) %>% 
  hc_exporting(enabled = TRUE) %>% 
  hc_legend(enabled = FALSE)


table.DownsideRisk(MVP_Comp_updated)

table.Drawdowns(MVP_Comp_updated)
SharpeRatio(MVP_Comp_updated)
table.Stats(MVP_Comp_updated)


############******  ***********



############****** Aufräumen Enviroment *************************************

# Liste aller Objekte im aktuellen Environment

all_objects <- ls()

# Liste der zu behaltenden Objekte
keep_objects <- c("all_mvp_returns_lagged_xts", "m_w_1", "m_w_2", "m_w_3", "m_w_4", "m_w_5",
                  "log_returns", "simple_returns", "dates", "MVP_Comp")

# Bestimmen der Objekte, die entfernt werden sollen
objects_to_remove <- setdiff(all_objects, keep_objects)

# Löschen der bestimmten Objekte
rm(list = objects_to_remove, envir = globalenv())

mvp_xts_1 <- na.omit(all_mvp_returns_lagged_xts[, "MVP_1"])
mvp_xts_2 <- na.omit(all_mvp_returns_lagged_xts[, "MVP_2"])
mvp_xts_3 <- na.omit(all_mvp_returns_lagged_xts[, "MVP_3"])
mvp_xts_4 <- na.omit(all_mvp_returns_lagged_xts[, "MVP_4"])
mvp_xts_5 <- na.omit(all_mvp_returns_lagged_xts[, "MVP_5"])

############****** DCC-GARCH *************************************

library(rmgarch)

fit_dcc_garch <- function(daily_returns) {
  spec <- ugarchspec(variance.model = list(model = "sGARCH"),
                     mean.model = list(armaOrder = c(1,1)),
                     distribution.model = "norm")
  dccspec <- dccspec(uspec = multispec(replicate(ncol(daily_returns), spec)),
                     dccOrder = c(1,1),
                     distribution = "mvnorm")
  dccfit <- dccfit(dccspec, data = daily_returns)
  return(dccfit)
}

# Funktion zur Berechnung der monatlichen MVP-Gewichte unter Verwendung des DCC-GARCH-Modells
calculate_monthly_mvp_weights_dcc <- function(log_returns, dcc_fit) {
  endpoints <- endpoints(log_returns, on = "months")
  monthly_weights <- list()
  
  for (i in 1:(length(endpoints) - 1)) {
    # Berechnen der dynamischen Kovarianzmatrix für den aktuellen Monat
    dcc_cov_matrix <- rcov(dcc_fit)[, , endpoints[i + 1]]
    
    # Berechnen der MVP-Gewichte
    mvp_weights <- minvar(dcc_cov_matrix, wmin = 0, wmax = 1)
    
    monthly_weights[[i]] <- mvp_weights
  }
  
  return(do.call(rbind, monthly_weights))
}

# Anpassen des DCC-GARCH-Modells auf tägliche Renditen
dcc_fit <- fit_dcc_garch(log_returns)

# Berechnen der monatlichen MVP-Gewichte mit DCC-GARCH
monthly_mvp_weights_dcc <- calculate_monthly_mvp_weights_dcc(log_returns, dcc_fit)


# Erstellen eines xts-Objekts mit den monatlichen Gewichten
monthly_mvp_weights_dcc_xts <- xts(monthly_mvp_weights_dcc, order.by = index(log_returns)[endpoints(log_returns, on = "months")[-1]])
colnames(monthly_mvp_weights_dcc_xts) <- asset_names


# Funktion zur Berechnung der monatlichen Renditen des MVP mit verzögerten Gewichten
calculate_monthly_mvp_returns_lagged <- function(weights_xts, log_returns_xts) {
  # Gewichte um einen Monat verzögern
  lagged_weights_xts <- lag.xts(monthly_mvp_weights_dcc_xts, k = +1)
  
  monthly_mvp_returns <- NULL
  
  # Gruppierung der Log-Renditen nach Monaten
  monthly_log_returns <- apply.monthly(log_returns, colSums)
  
  for (i in 2:nrow(lagged_weights_xts)) {
    # Extrahieren der verzögerten Gewichte
    lagged_weights <- as.numeric(lagged_weights_xts[i, ])
    
    # Extrahieren der monatlichen summierten Log-Renditen für den aktuellen Monat
    if (index(monthly_log_returns)[i] %in% index(log_returns)) {
      current_monthly_log_returns <- monthly_log_returns[index(lagged_weights_xts)[i], ]
      
      # Berechnen der monatlichen Rendite für das Portfolio
      monthly_return <- sum(current_monthly_log_returns * lagged_weights, na.rm = TRUE)
      monthly_mvp_returns <- c(monthly_mvp_returns, monthly_return)
    }
  }
  
  return(xts(monthly_mvp_returns, order.by = index(lagged_weights_xts)[-1]))
}


# Berechnen der monatlichen Renditen des MVP mit verzögerten Gewichten
monthly_mvp_returns_lagged_dcc <- calculate_monthly_mvp_returns_lagged(monthly_mvp_weights_dcc_xts, log_returns)
mvp_w_dcc <- monthly_mvp_weights_dcc_xts
# Aktualisieren des MVP-Vergleichsobjekts
MVP_Comp_updated1 <- merge.xts(MVP_Comp, monthly_mvp_returns_lagged_dcc)
colnames(MVP_Comp_updated1)[ncol(MVP_Comp_updated1)] <- "DCC-GARCH"

# Visualisieren der Portfolio-Leistung
charts.PerformanceSummary(MVP_Comp_updated1, colorset = 2:9)


chart.StackedBar(
  mvp_w_dcc,
  main = "kie",
  colorset = 2:9,
  space = 0.2,
  cex.axis = 0.8,
  cex.legend = 0.8,
  cex.lab = 1,
  cex.labels = 0.8,
  cex.main = 1,
  xaxis = TRUE,
  legend.loc = "under",
  element.color = "darkgray",
  unstacked = TRUE,
  xlab = "Date",
  ylab = "Value",
  ylim = NULL,
  date.format = "%b %y",
  major.ticks = "auto",
  minor.ticks = TRUE,
  las = 0,
  xaxis.labels = NULL,
)


par(mfrow = c(3, 2))






















