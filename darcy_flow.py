from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple


@dataclass(frozen=True)
class SeepageModelConfig:
    length: float = 60.0
    depth: float = 20.0
    nx: int = 81
    ny: int = 31
    upstream_head: float = 10.0
    downstream_head: float = 2.0
    hydraulic_conductivity: float = 1.0
    wall_x: float = 20.0
    wall_embedment: float = 12.0
    max_iterations: int = 10_000
    tolerance: float = 1e-4


def _make_grid(ny: int, nx: int, value: float = 0.0) -> List[List[float]]:
    return [[value for _ in range(nx)] for _ in range(ny)]


def solve_pressure_head(config: SeepageModelConfig) -> Tuple[Sequence[Sequence[float]], float]:
    if config.nx < 3 or config.ny < 3:
        raise ValueError("nx and ny must be at least 3")
    if config.length <= 0 or config.depth <= 0:
        raise ValueError("length and depth must be positive")
    if config.hydraulic_conductivity <= 0:
        raise ValueError("hydraulic_conductivity must be positive")

    dx = config.length / (config.nx - 1)
    dy = config.depth / (config.ny - 1)
    wall_i = max(1, min(config.nx - 2, round(config.wall_x / dx)))
    wall_j_max = max(0, min(config.ny - 1, round(config.wall_embedment / dy)))

    head = _make_grid(config.ny, config.nx, config.downstream_head)

    for j in range(config.ny):
        blend = 1.0 - (j / max(config.ny - 1, 1))
        initial = config.downstream_head + blend * (config.upstream_head - config.downstream_head)
        for i in range(config.nx):
            head[j][i] = initial
        head[j][0] = config.upstream_head
        head[j][-1] = config.downstream_head

    for _ in range(config.max_iterations):
        max_change = 0.0

        for j in range(config.ny):
            for i in range(1, config.nx - 1):
                if i == 0 or i == config.nx - 1:
                    continue
                if j == 0:
                    north = head[1][i]
                else:
                    north = head[j - 1][i]
                if j == config.ny - 1:
                    south = head[config.ny - 2][i]
                else:
                    south = head[j + 1][i]

                west = head[j][i - 1]
                east = head[j][i + 1]

                if j <= wall_j_max:
                    if i == wall_i:
                        west = head[j][i]
                    if i == wall_i - 1:
                        east = head[j][i]

                new_value = ((west + east) * dy * dy + (north + south) * dx * dx) / (2.0 * (dx * dx + dy * dy))
                change = abs(new_value - head[j][i])
                if change > max_change:
                    max_change = change
                head[j][i] = new_value

            head[j][0] = config.upstream_head
            head[j][-1] = config.downstream_head

        if max_change < config.tolerance:
            break

    outlet_gradient_sum = 0.0
    for j in range(config.ny):
        dh_dx = (head[j][-1] - head[j][-2]) / dx
        outlet_gradient_sum += dh_dx
    flow_rate = -config.hydraulic_conductivity * outlet_gradient_sum * dy

    return head, flow_rate


if __name__ == "__main__":
    cfg = SeepageModelConfig()
    _, outlet_q = solve_pressure_head(cfg)
    print(f"Estimated outlet seepage flow rate: {outlet_q:.6f}")
