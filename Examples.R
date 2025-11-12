### Use real RCTs with ordinal outcomes as pilot studies and compute the different sample sizes

##### Boulesteix, neonatology study, 2025
p_C_boulesteix <- c(0.265, 0.275, 0.247, 0.151, 0.02, 0.042)
p_E_boulesteix <- c(0.475, 0.18, 0.15, 0.137, 0.018, 0.04)

result_AfS_boulesteix <-samplesize_AfS(p_C=p_C_boulesteix, p_E=p_E_boulesteix, alpha=0.05, beta=0.2, r=1)
result_ttestord_boulesteix <-samplesize_ttestord(p_C=p_C_boulesteix, p_E=p_E_boulesteix, alpha=0.05, beta=0.2, r=1)
result_po_boulesteix <- samplesize_po_NN(p_C=p_C_boulesteix, p_E=p_E_boulesteix, alpha=0.05, beta=0.2, r=1)
data.frame("method"= c("WMW", "t-test", "PO"),
           "sample_size"=c(result_AfS_boulesteix$n_total, result_ttestord_boulesteix$n_total, result_po_boulesteix$n_total),
           "power" = c(result_AfS_boulesteix$actual_power, result_ttestord_boulesteix$actual_power, result_po_boulesteix$actual_power)
           )




###Baud O, Trousson C, Biran V, Leroy E, Mohamed D, Alberti C. Association
  # between early low-dose hydrocortisone therapy in extremely preterm
  # neonates and neurodevelopmental outcomes at 2 years of age. JAMA.
  # 2017;317(13):1329–37. https:// doi. org/ 10. 1001/ jama. 2017. 2692. United
  # States.
# (had to adjust the probabilities to sum to 1)
p_C_baud <- c(0.73, 0.2, 0.07)
p_E_baud <- c(0.70, 0.18, 0.12)
result_AfS_baud <-samplesize_AfS(p_C=p_C_baud, p_E=p_E_baud, alpha=0.05, beta=0.2, r=1)
result_ttestord_baud <-samplesize_ttestord(p_C=p_C_baud, p_E=p_E_baud, alpha=0.05, beta=0.2, r=1)
result_po_baud <- samplesize_po_NN(p_C=p_C_baud, p_E=p_E_baud, alpha=0.05, beta=0.2, r=1)
c("AfS"=result_AfS_baud$n_total, "ttest"=result_ttestord_baud$n_total, "po"=result_po_baud$n_total)
which.min(c("AfS"=result_AfS_baud$n_total, "ttest"=result_ttestord_baud$n_total, "po"=result_po_baud$n_total))
# t-test calculates much smaller sample size
#  man kann hald dann nicht die actual power berechnen, weil man die wahren prob. vektoren nicht kennt
# wie viel bringt das dann wirklich?
# oder soll ich die Wahrscheinlichkeiten als die Wahren nehmen?

### van den Berg LA, Dijkgraaf MG, Berkhemer OA, Fransen PS, Beumer D,
  # Lingsma HF, et al. Two-year outcome after endovascular treatment for
  # acute ischemic stroke. N Engl J Med. 2017;376(14):1341–9. https:// doi.
  # org/ 10. 1056/ NEJMo a1612 136. United States.
# after 2 years (had to adjust the probabilities to sum to 1)
p_C_berg <- c(0.026, 0.046, 0.299, 0.18, .062, 0.082, 0.305)
p_E_berg <- c(0.01, 0.051, 0.178, 0.168, 0.102, 0.107, 0.384)
result_AfS_berg <-samplesize_AfS(p_C=p_C_berg, p_E=p_E_berg, alpha=0.05, beta=0.2, r=1)
result_ttestord_berg <-samplesize_ttestord(p_C=p_C_berg, p_E=p_E_berg, alpha=0.05, beta=0.2, r=1)
result_po_berg <- samplesize_po_NN(p_C=p_C_berg, p_E=p_E_berg, alpha=0.05, beta=0.2, r=1)
c("AfS"=result_AfS_berg$n_total, "ttest"=result_ttestord_berg$n_total, "po"=result_po_berg$n_total)
which.min(c("AfS"=result_AfS_berg$n_total, "ttest"=result_ttestord_berg$n_total, "po"=result_po_berg$n_total))
# t-test a little smaller
