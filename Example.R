
# First run the file "Functions.R" to load the packages and function needed in this File

#### Use real RCTs with ordinal outcomes as pilot studies and compute the different sample sizes
#### Boulesteix (2025)
p_C_boulesteix <- c(0.265, 0.275, 0.247, 0.151, 0.02, 0.042)
p_E_boulesteix <- c(0.475, 0.18, 0.15, 0.137, 0.018, 0.04)

result_WMW_boulesteix <-samplesize_WMW(p_C=p_C_boulesteix, p_E=p_E_boulesteix, alpha=0.05, beta=0.2, r=1)
result_ttestord_boulesteix <-samplesize_ttestord(p_C=p_C_boulesteix, p_E=p_E_boulesteix, alpha=0.05, beta=0.2, r=1)
result_po_boulesteix <- samplesize_po_NN(p_C=p_C_boulesteix, p_E=p_E_boulesteix, alpha=0.05, beta=0.2, r=1)
data.frame("method"= c( "t-test", "WMW", "PO"),
           "sample_size"=c(result_ttestord_boulesteix$n_total, result_WMW_boulesteix$n_total, result_po_boulesteix$n_total),
           "power" = c(result_ttestord_boulesteix$actual_power, result_WMW_boulesteix$actual_power, result_po_boulesteix$actual_power)
)


