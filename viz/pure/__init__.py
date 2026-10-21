"""pure: a linguagem visual de madsondeluna, em Python.

Submodulos:

    palette   tokens de cor, sem dependencia externa
    mpl       tema matplotlib (requer matplotlib)
    plotly    template plotly (requer plotly)

Os submodulos de plotagem sao importados sob demanda, entao
`from pure import palette` funciona num ambiente sem matplotlib.
"""

from __future__ import annotations

from . import palette
from .palette import (
    CATEGORICAL,
    DIVERGING,
    ORDINAL_DARK,
    ORDINAL_LIGHT,
    SEQUENTIAL,
    ordinal,
    series,
    status,
    tokens,
)

__version__ = "1.6.0"

__all__ = [
    "palette",
    "CATEGORICAL",
    "DIVERGING",
    "SEQUENTIAL",
    "ORDINAL_LIGHT",
    "ORDINAL_DARK",
    "series",
    "tokens",
    "status",
    "ordinal",
    "__version__",
]


def __getattr__(name: str):
    if name in ("mpl", "plotly"):
        import importlib
        import sys

        module = importlib.import_module(f".{name}", __name__)
        setattr(sys.modules[__name__], name, module)
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
