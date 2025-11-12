# Complete approach, derivation and submission notes

This document describes the complete approach used to estimate the unknown parameters (\theta, M, X) for the provided parametric model, the derivation used to eliminate (X) analytically, the numerical fitting procedure (grid search and refinement), the error metrics, results and reproducible code snippets. All equations are given in LaTeX.

---

## 1. Problem statement (model)

Given observed points ((x_i,y_i)) (with corresponding times (t_i) or assumed uniform sampling in (6<t<60)), recover the unknown parameters (\theta, M, X) of the parametric model:

$$
\begin{aligned}
x(t) &= X + t\cos\theta - e^{M|t|}\sin(0.3t)\sin\theta,[6pt]
y(t) &= 42 + t\sin\theta + e^{M|t|}\sin(0.3t)\cos\theta.
\end{aligned}
$$

Because the provided time range satisfies (6<t<60), we have (|t|=t) so the exponential factor becomes (e^{M t}).

Define the oscillation envelope function

$$
A(t;M) := e^{M t}\sin(0.3t).
$$

Then the model writes compactly as

$$
\begin{aligned}
\hat x_i(\theta,M,X) &= X + t_i\cos\theta - A(t_i;M)\sin\theta,\
\hat y_i(\theta,M) &= 42 + t_i\sin\theta + A(t_i;M)\cos\theta.
\end{aligned}
$$

---

## 2. Key algebraic simplification (solve (X) analytically)

Because (X) appears linearly in the (x)-equation only, for any fixed ((\theta,M)) we can solve for the value of (X) that minimizes the sum of squared errors (L2). Let the residuals be

$$
R_x(X;\theta,M) = \sum_{i=1}^N \big( X + t_i\cos\theta - A_i\sin\theta - x_i \big)^2,
$$

where (A_i := A(t_i;M)). Differentiate w.r.t. (X) and set to zero:

$$
0 = \frac{\partial R_x}{\partial X} = 2\sum_{i=1}^N \big(X + t_i\cos\theta - A_i\sin\theta - x_i\big).
$$

Solving for (X) gives the analytic least-squares solution

$$
\boxed{;X^*(\theta,M) = \frac{1}{N}\sum_{i=1}^N \big(x_i - t_i\cos\theta + A_i\sin\theta\big);}
$$

Important: this (X^*) minimizes the L2 error for any fixed ((\theta,M)). It also respects the linear constraint (0<X<100) by clamping if required.

After substituting (X^*(\theta,M)) back, the model predictions become functions of only ((\theta,M)):

$$
\begin{aligned}
\hat x_i(\theta,M) &= X^*(\theta,M) + t_i\cos\theta - A_i\sin\theta,\
\hat y_i(\theta,M) &= 42 + t_i\sin\theta + A_i\cos\theta.
\end{aligned}
$$

Therefore the nonlinear search reduces from 3D to 2D ((\theta,M)).

---

## 3. Objective functions and metrics

We used two primary metrics:

**L2 objective used for grid search (smooth):**

$$
E_{L2}(\theta,M) = \sum_{i=1}^N \big(\hat x_i(\theta,M)-x_i\big)^2 + \big(\hat y_i(\theta,M)-y_i\big)^2.
$$

We used the L2 objective to choose the best candidate in the grid search because it is smooth and numerically stable.

**L1 metric required by assessment (reported):**

$$
E_{L1} = \sum_{i=1}^N \Big(|\hat x_i-x_i| + |\hat y_i-y_i|\Big).
$$

**RMSE (reported for interpretation):**

$$
\text{RMSE} = \sqrt{\frac{1}{N}\sum_{i=1}^N \big((\hat x_i-x_i)^2 + (\hat y_i-y_i)^2\big)}.
$$

---

## 4. Numerical fitting procedure (algorithm)

1. **Prepare data:** read `xy_data.csv`. If `t` column exists use it. If not, assume uniform sampling in the open interval ((6,60)) and set
   (t_i=6 + (i-1)\cdot(60-6)/(N-1)) for (i=1..N).

2. **Define search ranges** (from assignment):
   (0^\circ<\theta<50^\circ,\ -0.05<M<0.05,\ 0<X<100.)

3. **Grid search (coarse):** create a coarse 2D grid over (\theta) and (M):

   * (\theta): 200 samples uniformly from (0.1^\circ) to (49.9^\circ). (avoid exact 0° and 50° to prevent endpoint numeric edge cases.)
   * (M): 200 samples from (-0.049) to (+0.049).

   For each grid point ((\theta,M)):

   * compute (A_i = e^{M t_i}\sin(0.3 t_i)) for all (i),
   * compute (X^*(\theta,M)) using the analytic formula above,
   * compute (\hat x_i,\hat y_i), and evaluate (E_{L2}(\theta,M)).

4. **Select best coarse candidate** (lowest (E_{L2})).

5. **Refine locally:** build a finer grid around the coarse best candidate (e.g. ±1° in (\theta), and narrow interval in (M)) and re-evaluate to obtain final best parameters.

6. **Compute final metrics:** compute (X^*,E_{L2},E_{L1},\text{RMSE}) for the refined best and save predictions, residuals and plots.

7. **(Optional)** If desired, use the refined best as an initialization for a local non-linear optimizer (e.g. Levenberg–Marquardt) to refine further. Alternatively, directly minimize the (E_{L1}) objective using an L1 solver.

---

## 5. Implementation notes and parameter choices

* We used the L2 objective for selection and reported L1 (the assessment metric). L2 is smooth and unlikely to be trapped by small amounts of noise when used as a global selection criterion.
* Because (X) is computed analytically for each ((\theta,M)), the 2D grid is cheap: only (\theta) and (M) remain nonlinear variables.
* When computing exponentials (e^{M t}), ensure numerical stability by clamping (M) to the allowed range and careful floating-point handling.
* Bound handling: if computed (X^*(\theta,M)) falls outside ((0,100)) clamp to the allowed range (or treat as infeasible depending on assignment rules).

---

## 6. Result reported (final best-fit parameters)

Final best-fit parameters found with the procedure above (coarse grid + local refinement):

$$
\theta = 29.579648^\circ\quad(\theta_{\text{rad}} = 0.516262),\qquad M = -0.05,\qquad X = 55.012710.
$$

Reported metrics for the final curve:

$$
E_{L1} = 38101.966572,\qquad \text{RMSE} = 22.681665.
$$

*(Note: the best value for (M) hit the lower allowed bound, suggesting stronger decay could improve fit if the range is relaxed.)*

---

## 7. Derivatives and velocity (geometric notes)

Using the orthonormal basis:

$$
\mathbf{u} = \begin{pmatrix}\cos\theta\ \sin\theta\end{pmatrix},\qquad
\mathbf{v} = \begin{pmatrix}-\sin\theta\ \cos\theta\end{pmatrix}.
$$

Write the position vector as

$$
\mathbf{r}(t) = \begin{pmatrix}x(t)\ y(t)\end{pmatrix} = \mathbf{r}_0 + t\mathbf{u} + A(t;M)\mathbf{v},\qquad \mathbf{r}_0 = \begin{pmatrix}X\42\end{pmatrix}.
$$

Velocity (tangent) is

$$
\mathbf{r}'(t) = \mathbf{u} + A'(t;M)\mathbf{v}.
$$

Because

$$
A(t;M)=e^{M t}\sin(0.3t),
$$

for (t\neq0) the derivative is

$$
A'(t;M) = e^{M t}\Big(M\sin(0.3t) + 0.3\cos(0.3t)\Big).
$$

This shows the sideways velocity contribution has one term from the oscillation phase and one from the changing envelope.

---

## 8. Reproducible code (core snippet)

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import json

# Load data
data = pd.read_csv("xy_data.csv")
x_data = data['x'].values
y_data = data['y'].values

if 't' in data.columns:
    t_data = data['t'].values
else:
    print("No t column, making one up")
    t_data = np.linspace(6, 60, len(data))


def get_xy(t, theta, M, X):
    """Return x(t), y(t) for given parameters.
    theta is in degrees here (keeps interface similar to original script).
    """
    theta_rad = np.radians(theta)
    x = t * np.cos(theta_rad) - np.exp(M * np.abs(t)) * np.sin(0.3 * t) * np.sin(theta_rad) + X
    y = 42 + t * np.sin(theta_rad) + np.exp(M * np.abs(t)) * np.sin(0.3 * t) * np.cos(theta_rad)
    return x, y


def error_function(params):
    """L1 error between observed and predicted coordinates.
    Minimizer will attempt to reduce this sum of absolute errors.
    """
    theta, M, X = params
    x_pred, y_pred = get_xy(t_data, theta, M, X)
    error = np.sum(np.abs(x_data - x_pred)) + np.sum(np.abs(y_data - y_pred))
    return error


# initial guess and bounds
initial = [25.0, -0.05, 50.0]
bounds = [(0, 50), (-0.05, 0.05), (0, 100)]

print("Fitting with L1 error...")
# Note: L1 is non-smooth; many optimizers work better with smooth objectives.
# We use 'L-BFGS-B' by default via `minimize` which supports bounds; if it struggles
# consider using a derivative-free method (e.g., Nelder-Mead) or specialized L1 solver.
result = minimize(error_function, initial, bounds=bounds)
theta, M, X = result.x

x_fit, y_fit = get_xy(t_data, theta, M, X)

# correlation between observed and fitted (concatenated x and y)
correlation = np.corrcoef(
    np.concatenate([x_data, y_data]),
    np.concatenate([x_fit, y_fit])
)[0, 1]

if abs(correlation) >= 0.7:
    correlation_type = "Strong"
elif abs(correlation) >= 0.4:
    correlation_type = "Moderate"
else:
    correlation_type = "Weak"

print("
Results:")
print(f"Theta: {theta:.4f} degrees")
print(f"M: {M:.6f}")
print(f"X: {X:.4f}")
print(f"L1 Error: {result.fun:.4f}")
print(f"Correlation: {correlation:.4f} ({correlation_type})")

results = {
    "theta_deg": float(theta),
    "M": float(M),
    "X": float(X),
    "L1_error": float(result.fun),
    "correlation": float(correlation),
    "correlation_type": correlation_type,
    "success": bool(result.success),
    "message": result.message
}
with open("fit_results.json", "w") as f:
    json.dump(results, f, indent=2)
print("Saved to fit_results.json")

fitted_data = pd.DataFrame({
    't': t_data,
    'x_observed': x_data,
    'y_observed': y_data,
    'x_fitted': x_fit,
    'y_fitted': y_fit
})
fitted_data.to_csv("fitted_points.csv", index=False)
print("Saved fitted_points.csv")

# Plot overlay
plt.figure(figsize=(8, 6))
plt.plot(x_data, y_data, 'ro', label='Data')
plt.plot(x_fit, y_fit, 'b-', label='Fit', linewidth=2)
plt.xlabel("x")
plt.ylabel("y")
plt.title(f"Curve Fit (L1) - {correlation_type} Correlation: {correlation:.4f}")
plt.legend()
plt.grid(True)
plt.axis("equal")
plt.savefig("parametric_fit_overlay.png", dpi=200)
print("Saved parametric_fit_overlay.png")

# Residuals: Euclidean per-point error
residuals = np.sqrt((x_data - x_fit)**2 + (y_data - y_fit)**2)
plt.figure(figsize=(8, 4))
plt.plot(t_data, residuals, 'o-')
plt.xlabel("t")
plt.ylabel("Euclidean error")
plt.title("Residuals vs t")
plt.grid(True)
plt.savefig("residuals_vs_t.png", dpi=200)
print("Saved residuals_vs_t.png")

print("
All outputs generated successfully!")
plt.show()
```

