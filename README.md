# Research-and-Development_Flam
Fit a 2D parametric curve (θ, M, X) to provided (x,y) samples and report L1 error.

for the model
  x(t) = X + t cosθ − e^{M t} sin(0.3 t) sinθ
  y(t) = 42  + t sinθ + e^{M t} sin(0.3 t) cosθ
using the provided xy_data.csv (6 < t < 60). The repo contains the fitting code,
plots, fitted points and a short report showing the best-fit parameters and the
L1/RMSE metrics used for assessment.

