"""Tokens de cor da linguagem pure, em Python.

Espelha tokens/tokens.json. Nenhum valor e recalculado aqui: os hexes
sao os mesmos usados no CSS, para que uma figura exportada e a pagina
que a hospeda tenham exatamente a mesma tinta.

A paleta categorica passou os seis checks do validador de paletas nos
dois modos, com pior par adjacente delta-E 9.5 sob deuteranopia e 21.5
sob visao normal. A ordem dos slots e o mecanismo de seguranca: atribua
sempre em sequencia, nunca cicle.
"""

from __future__ import annotations

RAMP_SLATE: dict[str, str] = {
    "50": "#f4f6f9",
    "100": "#ebeef3",
    "150": "#e1e5ec",
    "200": "#dde1e9",
    "250": "#ced3dc",
    "300": "#adbcd0",
    "400": "#8ea2bc",
    "500": "#748cab",
    "600": "#667e9c",
    "700": "#3e5c76",
    "800": "#2a3f5c",
    "820": "#253752",
    "850": "#1d2d44",
    "900": "#141c2e",
    "950": "#0d1321",
}

RAMP_PAPER_LIKE: dict[str, str] = {
    "50": "#faf8f1",
    "100": "#f0ebd8",
    "150": "#e6e0ca",
    "250": "#d6cfb6",
}

RAMP_GRAPHITE: dict[str, str] = {
    "50": "#f3f6fb",
    "100": "#ebeef3",
    "150": "#e2e5ea",
    "200": "#dfe1e6",
    "250": "#d1d3d8",
    "300": "#b8babf",
    "400": "#9b9ea6",
    "500": "#868a92",
    "600": "#787c83",
    "700": "#555960",
    "800": "#3b3e45",
    "820": "#303339",
    "850": "#27282c",
    "900": "#18191d",
    "950": "#0e0f13",
}

LIGHT: dict[str, str] = {
    "bg": "#f4f6f9",
    "surface": "#ebeef3",
    "surface_hover": "#e1e5ec",
    "dim": "#dde1e9",
    "border": "#ced3dc",
    "border_hover": "#748cab",
    "text": "#0d1321",
    "muted": "#3e5c76",
    "accent": "#3e5c76",
    "secondary": "#748cab",
}

DEEP_BLUE: dict[str, str] = {
    "bg": "#0d1321",
    "surface": "#1d2d44",
    "surface_hover": "#253752",
    "dim": "#141c2e",
    "border": "#2a3f5c",
    "border_hover": "#3e5c76",
    "text": "#f4f6f9",
    # muted subiu para o passo 400 em 1.5.1: em #748cab ele ficava em
    # 4,02 sobre a superficie, abaixo do piso de 4,5 para texto pequeno,
    # e a camada de agente inteira poe texto de 12px em muted sobre cartao
    "muted": "#8ea2bc",
    # accent sobe junto, para o 300: precisa continuar um passo mais claro
    # que muted, como no modo dark (graphite-300 sobre graphite-400)
    "accent": "#adbcd0",
    "secondary": "#748cab",
}

PAPER_LIKE: dict[str, str] = {
    "bg": "#faf8f1",
    "surface": "#f0ebd8",
    "surface_hover": "#e6e0ca",
    "dim": "#e6e0ca",
    "border": "#d6cfb6",
    "border_hover": "#748cab",
    "text": "#0d1321",
    "muted": "#3e5c76",
    "accent": "#3e5c76",
    "secondary": "#748cab",
}

# Cinza chumbo a quase preto, sem azul. Todo par abre mais contraste que
# no deep-blue: muted sobe de 5,37 para 7,15 e secondary deixa de ser
# apenas decoracao, passando a 5,53.
DARK: dict[str, str] = {
    "bg": "#0e0f13",
    "surface": "#27282c",
    "surface_hover": "#303339",
    "dim": "#18191d",
    "border": "#3b3e45",
    "border_hover": "#555960",
    "text": "#f3f6fb",
    "muted": "#9b9ea6",
    "accent": "#b8babf",
    "secondary": "#868a92",
}

MODES: dict[str, dict[str, str]] = {
    "light": LIGHT,
    "paper-like": PAPER_LIKE,
    "deep-blue": DEEP_BLUE,
    "dark": DARK,
}

# Ordem fixa. O nono grupo de dados nunca ganha uma cor nova: vira "Outros",
# vira facetas, ou vira um segundo eixo de codificacao (marcador, textura).
CATEGORICAL: list[str] = [
    "#3973b1",  # 1 blue
    "#9f8322",  # 2 gold
    "#9e527f",  # 3 magenta
    "#4c985f",  # 4 green
    "#745ba5",  # 5 violet
    "#ba6f3e",  # 6 orange
    "#1990ad",  # 7 teal
    "#ac5551",  # 8 red
]

CATEGORICAL_NAMES: list[str] = [
    "blue", "gold", "magenta", "green", "violet", "orange", "teal", "red",
]

# Formas em que qualquer par de marcas pode encostar (dispersao, bolha,
# mapa, pequenos multiplos) so passam o teste de todos-os-pares com tres
# series. Acima disso: agrupe ou facete, nao troque a paleta.
ALL_PAIRS_SAFE_MAX = 3

# Rampa ordinal: passos discretos, uma so matiz, delta-L visivel.
# O escuro nao e o claro invertido: sao passos proprios, validados
# contra o fundo escuro.
ORDINAL_LIGHT: list[str] = [
    "#7bb1e9", "#5b9ddf", "#3b89d2", "#1f74bf", "#0f61a5", "#064f89", "#063d6b",
]
ORDINAL_DARK: list[str] = [
    "#b6d5f4", "#95c1ee", "#73ade7", "#5398dd", "#3284d0", "#1870ba", "#105d9d",
]

# Rampa sequencial continua, para heatmap e mapa de contato.
SEQUENTIAL: list[str] = [
    "#e5f0fc", "#cbe1f8", "#aacbee", "#82b3e5", "#5497d9",
    "#2b7ec9", "#1666aa", "#0e4f86", "#0e375c",
]

# Divergente: azul e ambar, cinza neutro no meio. Nunca uma matiz no centro.
DIVERGING: list[str] = [
    "#0e4f86", "#1f74bf", "#6da3da", "#b3ceeb", "#e2e5e8",
    "#eacab1", "#dc9a6c", "#c7692c", "#9e4421",
]

# Significado reservado. Nunca reaproveitar como "serie 5".
# Sempre acompanhado de icone ou rotulo, nunca cor sozinha.
STATUS_LIGHT: dict[str, str] = {
    "good": "#376e48",
    "warning": "#79601b",
    "serious": "#844e31",
    "critical": "#864544",
}
# os dois modos escuros compartilham as mesmas cores de status
STATUS_DARK: dict[str, str] = {
    "good": "#82bb90",
    "warning": "#c9ae6d",
    "serious": "#d89b7c",
    "critical": "#d18885",
}

FONT_SANS = ["Public Sans", "Helvetica Neue", "Arial", "DejaVu Sans"]
FONT_MONO = ["Spline Sans Mono", "Sometype Mono", "Menlo", "DejaVu Sans Mono"]


def series(n: int) -> list[str]:
    """Devolve as n primeiras cores categoricas, em ordem.

    Levanta ValueError acima de oito: a linguagem nao gera matiz nova.
    """
    if n < 1:
        raise ValueError("n deve ser >= 1")
    if n > len(CATEGORICAL):
        raise ValueError(
            f"{n} series pedidas, a paleta tem {len(CATEGORICAL)}. "
            "Agrupe o excedente em 'Outros', facete, ou use uma segunda "
            "codificacao (marcador, textura). Nao cicle as cores."
        )
    return CATEGORICAL[:n]


def tokens(mode: str = "light") -> dict[str, str]:
    """Tokens semanticos do modo pedido: light, paper-like, deep-blue ou dark."""
    if mode not in MODES:
        raise ValueError(f"modo desconhecido: {mode!r}. Use {list(MODES)}.")
    return dict(MODES[mode])


def status(mode: str = "light") -> dict[str, str]:
    """Cores de status do modo pedido. paper-like usa as do modo claro."""
    return dict(STATUS_DARK if mode in ("dark", "deep-blue") else STATUS_LIGHT)


def ordinal(n: int, mode: str = "light") -> list[str]:
    """n passos ordinais, amostrados da rampa do modo pedido."""
    ramp = ORDINAL_DARK if mode in ("dark", "deep-blue") else ORDINAL_LIGHT
    if n < 2:
        raise ValueError("uma rampa ordinal precisa de pelo menos 2 passos")
    if n > len(ramp):
        raise ValueError(f"{n} passos pedidos, a rampa tem {len(ramp)}")
    if n == len(ramp):
        return list(ramp)
    step = (len(ramp) - 1) / (n - 1)
    return [ramp[round(i * step)] for i in range(n)]
