# PIIM
Physics-Informed Neural Networks on Darcy Flow below a Water Reservoir

## Minimal seepage estimator

This repository now includes `darcy_flow.py`, a compact Darcy-flow pressure-head solver for flow through a homogeneous porous aquifer below a dam with:

- fixed upstream/downstream hydraulic heads,
- impermeable top/bottom boundaries (no-flow),
- an impermeable retaining wall with configurable embedment depth.

It computes:

1. the pressure-head field over a 2D aquifer grid, and
2. the estimated outlet seepage flow rate at the downstream boundary.

Run:

```bash
python darcy_flow.py
```
