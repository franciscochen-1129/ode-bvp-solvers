# ode-bvp-solvers

`ode-bvp-solvers` is a compact MATLAB repository for three standard numerical methods for two-point boundary value problems (BVPs):

- `finite_difference_bvp` for linear second-order BVPs
- `rayleigh_ritz_bvp` for variational / FEM-style solves
- `shooting_bvp_rk2` for shooting with RK2 time stepping and secant updates

The project is intentionally lightweight and classroom-friendly: each method is implemented in a small, readable MATLAB function, and the repository includes smoke tests and runnable examples.

## Repository Layout

```text
.
├── +odebvp/                 % package-style implementations
├── examples/                % minimal runnable demos
├── tests/                   % smoke checks
├── finite_difference_bvp.m  % standard public entry point
├── rayleigh_ritz_bvp.m      % standard public entry point
├── shooting_bvp_rk2.m       % standard public entry point
└── README.md
```

## Included Methods

### `finite_difference_bvp`

Solves linear BVPs of the form

```text
u''(x) + P(x)u'(x) + Q(x)u(x) = R(x),   a < x < b
u(a) = alpha,  u(b) = beta
```

using a second-order central-difference discretization on uniform grids.

### `rayleigh_ritz_bvp`

Solves

```text
-(P(x)v'(x))' + Q(x)v(x) = f(x),   a < x < b
v(a) = v(b) = 0
```

with piecewise linear basis functions on a uniform mesh.

### `shooting_bvp_rk2`

Solves nonlinear or linear second-order scalar BVPs

```text
y'' = F(x, y, y')
y(a) = alpha,  y(b) = beta
```

using:

- shooting on the unknown initial slope
- RK2 / improved Euler for the IVP solve
- secant iteration for boundary matching

## Quick Start

Clone `ode-bvp-solvers` and add it to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

### Finite Difference Example

```matlab
P = @(x) 0;
Q = @(x) -1;
R = @(x) -x;
[N_values, h_values, sols] = finite_difference_bvp(P, Q, R, 0, 1, 0, 0, 8, 16);
disp([N_values(:), h_values(:)]);
plot(sols{end}.x, sols{end}.u, '-o');
```

### Rayleigh-Ritz Example

```matlab
P = @(x) 1 + 0*x;
Q = @(x) 0*x;
f = @(x) pi^2 * sin(pi*x);
[N_values, h_values, sols] = rayleigh_ritz_bvp(P, Q, f, 0, 1, 8, 16);
plot(sols{end}.x, sols{end}.v, '-o');
```

### Shooting Example

```matlab
F = @(x, y, yp) -y;
[x, y, s_final, it] = shooting_bvp_rk2(F, 0, pi/2, 0, 1, 100, 0.5, 1.5, 1e-8, 25);
plot(x, y, 'LineWidth', 1.5);
```

## Public API

Main entry points:

- `finite_difference_bvp(...)`
- `rayleigh_ritz_bvp(...)`
- `shooting_bvp_rk2(...)`

These top-level functions delegate to the package implementations in `+odebvp/`.

## Examples And Tests

Run the demos:

```matlab
run('examples/demo_finite_difference.m');
run('examples/demo_rayleigh_ritz.m');
run('examples/demo_shooting.m');
```

Run the smoke tests:

```matlab
run('tests/smoke_test.m');
```

## Notes

- All methods use uniform meshes.
- The code favors clarity over heavy optimization.
- No external MATLAB toolboxes are required beyond standard numerical routines such as `integral`.
