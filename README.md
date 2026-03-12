# Stub---Casamento-de-impedância
Este projeto consiste em um programa em Python para análise e cálculo de casamento de impedância por stub em linhas de transmissão. O objetivo é automatizar os cálculos mais comuns da disciplina de Ondas e Linhas, permitindo que o usuário informe apenas os parâmetros do circuito e obtenha as soluções de forma rápida e organizada.


"""
Calculadora de casamento por stub para linhas de transmissão.

Casos implementados:
- Single-stub shunt (derivação) com stub em curto
- Single-stub shunt (derivação) com stub em aberto
- Single-stub series (série) com stub em curto
- Single-stub series (série) com stub em aberto
- Double-stub shunt (derivação) com stubs em curto/aberto
- Double-stub series (série) com stubs em curto/aberto

Entradas:
- Z0: impedância característica da linha [ohms]
- ZL: impedância da carga [ohms] (pode ser complexa)
- f: frequência [Hz]
- vp: velocidade de propagação [m/s] (padrão: 3e8)
- arrangement: "shunt" ou "series"
- stub_termination: "short" ou "open"
- stub_spacing_m: espaçamento entre os dois stubs [m] para double-stub

Saídas:
- soluções possíveis para:
    * distância da carga até o stub (single)
    * comprimento do stub (single)
    * comprimentos dos dois stubs (double)
    * impedância/admitância normalizada
    * coeficiente de reflexão e VSWR da carga
"""

import math
from dataclasses import dataclass
from typing import List, Optional


# ============================================================
# Utilidades
# ============================================================

def cfmt(z: complex, ndigits: int = 4) -> str:
    a = round(z.real, ndigits)
    b = round(z.imag, ndigits)
    if abs(b) < 10**(-ndigits):
        return f"{a}"
    sign = "+" if b >= 0 else "-"
    return f"{a} {sign} j{abs(b)}"


def gamma_from_impedance(ZL: complex, Z0: float) -> complex:
    return (ZL - Z0) / (ZL + Z0)


def vswr_from_gamma(gamma: complex) -> float:
    mag = abs(gamma)
    if math.isclose(mag, 1.0, rel_tol=1e-12, abs_tol=1e-12):
        return math.inf
    return (1 + mag) / (1 - mag)


def wavelength(f: float, vp: float = 3e8) -> float:
    return vp / f


def beta_from_f(f: float, vp: float = 3e8) -> float:
    return 2 * math.pi / wavelength(f, vp)


def normalize_impedance(ZL: complex, Z0: float) -> complex:
    return ZL / Z0


def normalize_admittance(ZL: complex, Z0: float) -> complex:
    return 1 / normalize_impedance(ZL, Z0)


def wrap_length_0_to_half_lambda(length: float, lam: float) -> float:
    half = lam / 2
    length = length % half
    if length < 0:
        length += half
    return length


def line_input_impedance_normalized(zL: complex, beta: float, d: float) -> complex:
    t = math.tan(beta * d)
    return (zL + 1j * t) / (1 + 1j * zL * t)


def line_input_admittance_normalized(yL: complex, beta: float, d: float) -> complex:
    t = math.tan(beta * d)
    return (yL + 1j * t) / (1 + 1j * yL * t)


def shunt_stub_susceptance_from_length(beta: float, l: float, stub_termination: str) -> float:
    if stub_termination == "short":
        # y_stub = -j cot(beta l)
        s = math.sin(beta * l)
        c = math.cos(beta * l)
        if abs(s) < 1e-12:
            return math.copysign(math.inf, -c if c != 0 else 1.0)
        return -(c / s)
    elif stub_termination == "open":
        return math.tan(beta * l)
    raise ValueError("stub_termination deve ser 'short' ou 'open'.")


def series_stub_reactance_from_length(beta: float, l: float, stub_termination: str) -> float:
    if stub_termination == "short":
        return math.tan(beta * l)
    elif stub_termination == "open":
        s = math.sin(beta * l)
        c = math.cos(beta * l)
        if abs(s) < 1e-12:
            return math.copysign(math.inf, -c if c != 0 else 1.0)
        return -(c / s)
    raise ValueError("stub_termination deve ser 'short' ou 'open'.")


def shunt_length_from_susceptance(b_needed: float, beta: float, lam: float, stub_termination: str) -> float:
    if stub_termination == "short":
        if math.isclose(b_needed, 0.0, abs_tol=1e-14):
            l = lam / 4
        else:
            l = math.atan(-1 / b_needed) / beta
    else:
        l = math.atan(b_needed) / beta
    return wrap_length_0_to_half_lambda(l, lam)


def series_length_from_reactance(x_needed: float, beta: float, lam: float, stub_termination: str) -> float:
    if stub_termination == "short":
        l = math.atan(x_needed) / beta
    else:
        if math.isclose(x_needed, 0.0, abs_tol=1e-14):
            l = lam / 4
        else:
            l = math.atan(-1 / x_needed) / beta
    return wrap_length_0_to_half_lambda(l, lam)


# ============================================================
# Estruturas
# ============================================================

@dataclass
class StubSolution:
    arrangement: str
    stub_termination: str
    distance_from_load_m: float
    distance_from_load_lambda: float
    stub_length_m: float
    stub_length_lambda: float
    transformed_normalized_value: complex
    stub_normalized_reactive_value: float


@dataclass
class DoubleStubSolution:
    arrangement: str
    stub_termination: str
    first_stub_length_m: float
    first_stub_length_lambda: float
    second_stub_length_m: float
    second_stub_length_lambda: float
    stub_spacing_m: float
    stub_spacing_lambda: float
    first_stub_normalized_reactive_value: float
    second_stub_normalized_reactive_value: float
    normalized_value_after_first_stub: complex
    normalized_value_before_second_stub: complex


# ============================================================
# Single-stub
# ============================================================

def solve_single_stub_shunt(Z0: float, ZL: complex, f: float, vp: float = 3e8,
                            stub_termination: str = "short") -> List[StubSolution]:
    beta = beta_from_f(f, vp)
    lam = wavelength(f, vp)

    zL = normalize_impedance(ZL, Z0)
    yL = 1 / zL
    gL = yL.real
    bL = yL.imag

    A = (gL**2 + bL**2 - gL)
    B = -2 * bL
    C = (1 - gL)

    if math.isclose(A, 0.0, abs_tol=1e-14):
        ts = [] if math.isclose(B, 0.0, abs_tol=1e-14) else [-C / B]
    else:
        disc = B**2 - 4 * A * C
        if disc < -1e-12:
            ts = []
        else:
            disc = max(disc, 0.0)
            sqrt_disc = math.sqrt(disc)
            ts = [(-B + sqrt_disc) / (2 * A), (-B - sqrt_disc) / (2 * A)]

    solutions = []
    used = []

    for t in ts:
        d = math.atan(t) / beta
        if d < 0:
            d += lam / 2
        if any(abs(d - u) < 1e-10 for u in used):
            continue
        used.append(d)

        y_d = line_input_admittance_normalized(yL, beta, d)
        b_needed = -y_d.imag
        l = shunt_length_from_susceptance(b_needed, beta, lam, stub_termination)

        solutions.append(StubSolution(
            arrangement="shunt",
            stub_termination=stub_termination,
            distance_from_load_m=d,
            distance_from_load_lambda=d / lam,
            stub_length_m=l,
            stub_length_lambda=l / lam,
            transformed_normalized_value=y_d,
            stub_normalized_reactive_value=b_needed,
        ))

    solutions.sort(key=lambda s: s.distance_from_load_m)
    return solutions


def solve_single_stub_series(Z0: float, ZL: complex, f: float, vp: float = 3e8,
                             stub_termination: str = "short") -> List[StubSolution]:
    beta = beta_from_f(f, vp)
    lam = wavelength(f, vp)

    zL = normalize_impedance(ZL, Z0)
    rL = zL.real
    xL = zL.imag

    A = (rL - 1)
    B = -2 * xL
    C = (rL**2 + xL**2 - rL)

    if math.isclose(A, 0.0, abs_tol=1e-14):
        ts = [] if math.isclose(B, 0.0, abs_tol=1e-14) else [-C / B]
    else:
        disc = B**2 - 4 * A * C
        if disc < -1e-12:
            ts = []
        else:
            disc = max(disc, 0.0)
            sqrt_disc = math.sqrt(disc)
            ts = [(-B + sqrt_disc) / (2 * A), (-B - sqrt_disc) / (2 * A)]

    solutions = []
    used = []

    for t in ts:
        d = math.atan(t) / beta
        if d < 0:
            d += lam / 2
        if any(abs(d - u) < 1e-10 for u in used):
            continue
        used.append(d)

        z_d = line_input_impedance_normalized(zL, beta, d)
        x_needed = -z_d.imag
        l = series_length_from_reactance(x_needed, beta, lam, stub_termination)

        solutions.append(StubSolution(
            arrangement="series",
            stub_termination=stub_termination,
            distance_from_load_m=d,
            distance_from_load_lambda=d / lam,
            stub_length_m=l,
            stub_length_lambda=l / lam,
            transformed_normalized_value=z_d,
            stub_normalized_reactive_value=x_needed,
        ))

    solutions.sort(key=lambda s: s.distance_from_load_m)
    return solutions


# ============================================================
# Double-stub shunt
# ============================================================

def solve_double_stub_shunt(Z0: float, ZL: complex, f: float, stub_spacing_m: float,
                            vp: float = 3e8, stub_termination: str = "short") -> List[DoubleStubSolution]:
    beta = beta_from_f(f, vp)
    lam = wavelength(f, vp)
    d = stub_spacing_m

    zL = normalize_impedance(ZL, Z0)
    yL = 1 / zL
    g = yL.real
    b = yL.imag
    t = math.tan(beta * d)

    # Seja B1 a susceptância do primeiro stub.
    # y1 = g + j(b + B1)
    # Após propagar d até o 2º stub: y2 = (y1 + jt)/(1 + jy1 t)
    # Impor Re{y2}=1 gera:
    # a*u^2 + bq*u + c = 0, com u = (b + B1)
    # a = t^2
    # bq = -2t
    # c = g^2(1+t^2) - g
    a = t**2
    bq = -2 * t
    c = g**2 * (1 + t**2) - g

    if math.isclose(a, 0.0, abs_tol=1e-14):
        # Espaçamento ~ n*lambda/2 degenera o double-stub
        return []

    disc = bq**2 - 4 * a * c
    if disc < -1e-12:
        return []

    disc = max(disc, 0.0)
    sqrt_disc = math.sqrt(disc)
    u_roots = [(-bq + sqrt_disc) / (2 * a), (-bq - sqrt_disc) / (2 * a)]

    sols = []
    used = []

    for u in u_roots:
        B1 = u - b
        y1 = g + 1j * (b + B1)
        y2 = line_input_admittance_normalized(y1, beta, d)
        B2 = -y2.imag

        l1 = shunt_length_from_susceptance(B1, beta, lam, stub_termination)
        l2 = shunt_length_from_susceptance(B2, beta, lam, stub_termination)

        key = (round(l1, 12), round(l2, 12))
        if key in used:
            continue
        used.append(key)

        sols.append(DoubleStubSolution(
            arrangement="shunt",
            stub_termination=stub_termination,
            first_stub_length_m=l1,
            first_stub_length_lambda=l1 / lam,
            second_stub_length_m=l2,
            second_stub_length_lambda=l2 / lam,
            stub_spacing_m=d,
            stub_spacing_lambda=d / lam,
            first_stub_normalized_reactive_value=B1,
            second_stub_normalized_reactive_value=B2,
            normalized_value_after_first_stub=y1,
            normalized_value_before_second_stub=y2,
        ))

    sols.sort(key=lambda s: (s.first_stub_length_m, s.second_stub_length_m))
    return sols


# ============================================================
# Double-stub series
# ============================================================

def solve_double_stub_series(Z0: float, ZL: complex, f: float, stub_spacing_m: float,
                             vp: float = 3e8, stub_termination: str = "short") -> List[DoubleStubSolution]:
    beta = beta_from_f(f, vp)
    lam = wavelength(f, vp)
    d = stub_spacing_m

    zL = normalize_impedance(ZL, Z0)
    r = zL.real
    x = zL.imag
    t = math.tan(beta * d)

    # Seja X1 a reatância do primeiro stub.
    # z1 = r + j(x + X1)
    # Após propagar d: z2 = (z1 + jt)/(1 + j z1 t)
    # Impor Re{z2}=1 resulta em:
    # a*u^2 + bq*u + c = 0, com u = x + X1
    # a = 1 - r
    # bq = -2t
    # c = r - r^2 - r^2 t^2
    a = (1 - r)
    bq = -2 * t
    c = r - r**2 - r**2 * t**2

    if math.isclose(a, 0.0, abs_tol=1e-14):
        if math.isclose(bq, 0.0, abs_tol=1e-14):
            return []
        u_roots = [-c / bq]
    else:
        disc = bq**2 - 4 * a * c
        if disc < -1e-12:
            return []
        disc = max(disc, 0.0)
        sqrt_disc = math.sqrt(disc)
        u_roots = [(-bq + sqrt_disc) / (2 * a), (-bq - sqrt_disc) / (2 * a)]

    sols = []
    used = []

    for u in u_roots:
        X1 = u - x
        z1 = r + 1j * (x + X1)
        z2 = line_input_impedance_normalized(z1, beta, d)
        X2 = -z2.imag

        l1 = series_length_from_reactance(X1, beta, lam, stub_termination)
        l2 = series_length_from_reactance(X2, beta, lam, stub_termination)

        key = (round(l1, 12), round(l2, 12))
        if key in used:
            continue
        used.append(key)

        sols.append(DoubleStubSolution(
            arrangement="series",
            stub_termination=stub_termination,
            first_stub_length_m=l1,
            first_stub_length_lambda=l1 / lam,
            second_stub_length_m=l2,
            second_stub_length_lambda=l2 / lam,
            stub_spacing_m=d,
            stub_spacing_lambda=d / lam,
            first_stub_normalized_reactive_value=X1,
            second_stub_normalized_reactive_value=X2,
            normalized_value_after_first_stub=z1,
            normalized_value_before_second_stub=z2,
        ))

    sols.sort(key=lambda s: (s.first_stub_length_m, s.second_stub_length_m))
    return sols


# ============================================================
# Interfaces principais
# ============================================================

def solve_stub_matching(Z0: float, ZL: complex, f: float, vp: float = 3e8,
                        arrangement: str = "shunt",
                        stub_termination: str = "short") -> dict:
    arrangement = arrangement.lower().strip()
    stub_termination = stub_termination.lower().strip()

    gamma = gamma_from_impedance(ZL, Z0)
    vswr = vswr_from_gamma(gamma)
    lam = wavelength(f, vp)
    beta = beta_from_f(f, vp)

    if arrangement == "shunt":
        sols = solve_single_stub_shunt(Z0, ZL, f, vp, stub_termination)
    elif arrangement == "series":
        sols = solve_single_stub_series(Z0, ZL, f, vp, stub_termination)
    else:
        raise ValueError("arrangement deve ser 'shunt' ou 'series'.")

    return {
        "mode": "single",
        "Z0_ohms": Z0,
        "ZL_ohms": ZL,
        "frequency_hz": f,
        "vp_m_s": vp,
        "wavelength_m": lam,
        "beta_rad_m": beta,
        "gamma_load": gamma,
        "abs_gamma": abs(gamma),
        "vswr": vswr,
        "arrangement": arrangement,
        "stub_termination": stub_termination,
        "solutions": sols,
    }


def solve_double_stub_matching(Z0: float, ZL: complex, f: float, stub_spacing_m: float,
                               vp: float = 3e8, arrangement: str = "shunt",
                               stub_termination: str = "short") -> dict:
    arrangement = arrangement.lower().strip()
    stub_termination = stub_termination.lower().strip()

    gamma = gamma_from_impedance(ZL, Z0)
    vswr = vswr_from_gamma(gamma)
    lam = wavelength(f, vp)
    beta = beta_from_f(f, vp)

    if arrangement == "shunt":
        sols = solve_double_stub_shunt(Z0, ZL, f, stub_spacing_m, vp, stub_termination)
    elif arrangement == "series":
        sols = solve_double_stub_series(Z0, ZL, f, stub_spacing_m, vp, stub_termination)
    else:
        raise ValueError("arrangement deve ser 'shunt' ou 'series'.")

    return {
        "mode": "double",
        "Z0_ohms": Z0,
        "ZL_ohms": ZL,
        "frequency_hz": f,
        "vp_m_s": vp,
        "wavelength_m": lam,
        "beta_rad_m": beta,
        "gamma_load": gamma,
        "abs_gamma": abs(gamma),
        "vswr": vswr,
        "arrangement": arrangement,
        "stub_termination": stub_termination,
        "stub_spacing_m": stub_spacing_m,
        "stub_spacing_lambda": stub_spacing_m / lam,
        "solutions": sols,
    }


def print_report(result: dict) -> None:
    print("=" * 78)
    print(f"RELATÓRIO DE CASAMENTO POR STUB ({result['mode'].upper()})")
    print("=" * 78)
    print(f"Z0 = {result['Z0_ohms']} ohms")
    print(f"ZL = {cfmt(result['ZL_ohms'])} ohms")
    print(f"f  = {result['frequency_hz']} Hz")
    print(f"vp = {result['vp_m_s']} m/s")
    print(f"lambda = {result['wavelength_m']:.6f} m")
    print(f"beta   = {result['beta_rad_m']:.6f} rad/m")
    print(f"Gamma_L = {cfmt(result['gamma_load'])}")
    print(f"|Gamma_L| = {result['abs_gamma']:.6f}")
    print(f"VSWR = {result['vswr']:.6f}")
    print(f"Configuração: {result['arrangement']} / {result['stub_termination']}")
    if result["mode"] == "double":
        print(f"Espaçamento entre stubs = {result['stub_spacing_m']:.6f} m "
              f"= {result['stub_spacing_lambda']:.6f} λ")
    print("-" * 78)

    sols = result["solutions"]
    if not sols:
        print("Nenhuma solução encontrada para os parâmetros informados.")
        print("No caso de double-stub, isso pode ocorrer por restrição geométrica")
        print("do espaçamento entre stubs para essa carga específica.")
        print("-" * 78)
        return

    for i, s in enumerate(sols, start=1):
        print(f"Solução {i}:")
        if result["mode"] == "single":
            if s.arrangement == "shunt":
                nome = "Admitância normalizada no ponto"
                reativo = "b_stub necessário"
            else:
                nome = "Impedância normalizada no ponto"
                reativo = "x_stub necessário"

            print(f"  Distância da carga até o stub = {s.distance_from_load_m:.6f} m "
                  f"= {s.distance_from_load_lambda:.6f} λ")
            print(f"  Comprimento do stub          = {s.stub_length_m:.6f} m "
                  f"= {s.stub_length_lambda:.6f} λ")
            print(f"  {nome} = {cfmt(s.transformed_normalized_value)}")
            print(f"  {reativo} = {s.stub_normalized_reactive_value:.6f}")
        else:
            if s.arrangement == "shunt":
                nome1 = "Admitância após 1º stub"
                nome2 = "Admitância antes do 2º stub"
                re1 = "b1"
                re2 = "b2"
            else:
                nome1 = "Impedância após 1º stub"
                nome2 = "Impedância antes do 2º stub"
                re1 = "x1"
                re2 = "x2"

            print(f"  Comprimento do 1º stub = {s.first_stub_length_m:.6f} m "
                  f"= {s.first_stub_length_lambda:.6f} λ")
            print(f"  Comprimento do 2º stub = {s.second_stub_length_m:.6f} m "
                  f"= {s.second_stub_length_lambda:.6f} λ")
            print(f"  Espaçamento entre stubs = {s.stub_spacing_m:.6f} m "
                  f"= {s.stub_spacing_lambda:.6f} λ")
            print(f"  {re1} = {s.first_stub_normalized_reactive_value:.6f}")
            print(f"  {re2} = {s.second_stub_normalized_reactive_value:.6f}")
            print(f"  {nome1} = {cfmt(s.normalized_value_after_first_stub)}")
            print(f"  {nome2} = {cfmt(s.normalized_value_before_second_stub)}")
        print("-" * 78)


# ============================================================
# Exemplos
# ============================================================

if __name__ == "__main__":
    Z0 = 50
    ZL = 30 + 1j * 40
    f = 1e9
    vp = 3e8

    print("\nEXEMPLO 1 -> single stub shunt curto\n")
    res1 = solve_stub_matching(Z0, ZL, f, vp, arrangement="shunt", stub_termination="short")
    print_report(res1)

    print("\nEXEMPLO 2 -> single stub series aberto\n")
    res2 = solve_stub_matching(Z0, ZL, f, vp, arrangement="series", stub_termination="open")
    print_report(res2)

    print("\nEXEMPLO 3 -> double stub shunt curto, espaçamento = lambda/8\n")
    res3 = solve_double_stub_matching(
        Z0, ZL, f, wavelength(f, vp) / 8, vp,
        arrangement="shunt", stub_termination="short"
    )
    print_report(res3)

    print("\nEXEMPLO 4 -> double stub series aberto, espaçamento = lambda/8\n")
    res4 = solve_double_stub_matching(
        Z0, ZL, f, wavelength(f, vp) / 8, vp,
        arrangement="series", stub_termination="open"
    )
    print_report(res4)
