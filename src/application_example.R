# Application example
# This script applies the proposed methods and an existing method to a real dataset and compares the results.

# Required packages
# install.packages("faraway")
# install.packages("mada")
# install.packages("patchwork")
# install.packages("ggplot2")

source("src/functions/robustDTAMA.R")

path_dat = "data"
out_dir = "output/application_example"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

set.seed(1234)

# MMSE data (Arevalo et al.)
## 8 studies for diagnosis of conversion from MCI to ADD
dat = read.csv(file.path(path_dat,"MMSE.csv"))

# existing method (Bivariate Normal-Normal model with REML)
res_bnn = fit_bnn(Data = dat, method = "reml") 

# proposed method with alpha selected by the data-adaptive criterion based on the Hyvarinen score (alpha_H)
alpha_u_lim = 0.5
res_bnnDPD_alpha_H = fit_bnnDPD_grid_search(dat, gamma=0.01, alpha_lim=c(0.01, alpha_u_lim), verbose=TRUE)
# proposed method with fixed alpha = 0.33 (alpha_GES)
alpha_const = 0.33
res_bnnDPD_alpha_GES = fit_bnnDPD(dat, alpha=alpha_const, verbose=FALSE)

# Display Se, Sp, CIs, HKSJ CIs, and alpha for selected methods
extract_SeSp_CI <- function(method_name, se_sp, cis, hksj_cis = c(NA, NA, NA, NA), alpha = NA) {
  data.frame(
    method = method_name,
    Se = se_sp[1],
    Se_ci_l = cis[1],
    Se_ci_u = cis[2],
    Se_hksj_ci_l = hksj_cis[1],
    Se_hksj_ci_u = hksj_cis[2],
    Sp = se_sp[2],
    Sp_ci_l = cis[3],
    Sp_ci_u = cis[4],
    Sp_hksj_ci_l = hksj_cis[3],
    Sp_hksj_ci_u = hksj_cis[4],
    alpha = alpha,
    row.names = NULL
  )
}

res_SeSp_CI <- rbind(
  extract_SeSp_CI("BNN", res_bnn$SeSp, res_bnn$CIs, alpha = NA),
  extract_SeSp_CI(
    "BNNDPD (alpha_H)",
    res_bnnDPD_alpha_H$SeSp,
    c(res_bnnDPD_alpha_H$CIs[1, "ci_lower"],
      res_bnnDPD_alpha_H$CIs[1, "ci_upper"],
      res_bnnDPD_alpha_H$CIs[2, "ci_lower"],
      res_bnnDPD_alpha_H$CIs[2, "ci_upper"]),
    hksj_cis = c(
      res_bnnDPD_alpha_H$est_summary$ci_lower_hksj[6],
      res_bnnDPD_alpha_H$est_summary$ci_upper_hksj[6],
      res_bnnDPD_alpha_H$est_summary$ci_lower_hksj[7],
      res_bnnDPD_alpha_H$est_summary$ci_upper_hksj[7]
    ),
    alpha = res_bnnDPD_alpha_H$best_alpha
  ),
  extract_SeSp_CI(
    "BNNDPD (alpha_GES)",
    res_bnnDPD_alpha_GES$SeSp,
    c(res_bnnDPD_alpha_GES$CIs[1, "ci_lower"],
      res_bnnDPD_alpha_GES$CIs[1, "ci_upper"],
      res_bnnDPD_alpha_GES$CIs[2, "ci_lower"],
      res_bnnDPD_alpha_GES$CIs[2, "ci_upper"]),
    hksj_cis = c(
      res_bnnDPD_alpha_GES$est_summary$ci_lower_hksj[6],
      res_bnnDPD_alpha_GES$est_summary$ci_upper_hksj[6],
      res_bnnDPD_alpha_GES$est_summary$ci_lower_hksj[7],
      res_bnnDPD_alpha_GES$est_summary$ci_upper_hksj[7]
    ),
    alpha = alpha_const
  )
)

num_cols <- setdiff(names(res_SeSp_CI), "method")
res_SeSp_CI[, num_cols] <- round(res_SeSp_CI[, num_cols], 3)
print(res_SeSp_CI)

# plot of study-specific contribution measures
p=plot_dpd_contribution_rate(
  res_obj = res_bnnDPD_alpha_GES,
  ylim_g = c(0, 1),
  ylim_u = c(0, 0.4)
)
ggplot2::ggsave(
  filename = file.path(out_dir, "DPD_contribution.tiff"),
  plot = p$p_panel, width = 18, height = 14, units = "cm", dpi = 300
)

cat("All analyses completed!\n")
