# Stub - Casamento de Impedância

Este projeto consiste em um programa em **Python** para análise e cálculo de **casamento de impedância por stub** em **linhas de transmissão**. O objetivo é automatizar os cálculos mais comuns da disciplina de **Ondas e Linhas**, permitindo que o usuário informe apenas os parâmetros do circuito e obtenha as soluções de forma rápida e organizada.

## Funcionalidades

O programa resolve os seguintes casos:

- **Single-stub shunt** (derivação) com stub em curto
- **Single-stub shunt** (derivação) com stub em aberto
- **Single-stub series** (série) com stub em curto
- **Single-stub series** (série) com stub em aberto
- **Double-stub shunt** (derivação) com stubs em curto ou aberto
- **Double-stub series** (série) com stubs em curto ou aberto

## Parâmetros de entrada

O usuário informa:

- `Z0`: impedância característica da linha \[ohms]
- `ZL`: impedância da carga \[ohms] (pode ser complexa)
- `f`: frequência \[Hz]
- `vp`: velocidade de propagação \[m/s]
- `arrangement`: `"shunt"` ou `"series"`
- `stub_termination`: `"short"` ou `"open"`
- `stub_spacing_m`: espaçamento entre os dois stubs \[m] para o caso de **double-stub**

## Resultados fornecidos

O programa calcula e exibe:

- comprimento de onda \( \lambda \)
- constante de fase \( \beta \)
- coeficiente de reflexão da carga \( \Gamma \)
- módulo de \( \Gamma \)
- **VSWR**
- distância da carga até o stub no caso de **single-stub**
- comprimento do stub no caso de **single-stub**
- comprimentos do primeiro e do segundo stub no caso de **double-stub**
- impedância ou admitância normalizada em pontos importantes da linha
- soluções possíveis para o casamento

## Observações

- No caso de **double-stub**, o espaçamento entre os stubs deve ser fornecido.
- Dependendo da carga e do espaçamento, pode acontecer de **não existir solução**.
- Os comprimentos retornados ficam no intervalo:

```text
0 <= l < λ/2
