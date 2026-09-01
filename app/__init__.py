"""The Fermentation Console — the interface layer.

Sits *above* ``scenario``/``validation`` in the dependency order: it imports the engine and
the engine never imports it. Its dependencies (Streamlit, Plotly) live in the ``ui``
dependency group, so ``uv sync`` alone still installs a four-package research library.

``render.py`` is deliberately framework-free — it turns a finished run into figures and
panels and knows nothing about Streamlit — which is what lets the same code back both the
live app (``main.py``) and the standalone written report (``report.py``).
"""
