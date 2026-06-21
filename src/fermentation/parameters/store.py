"""Loading and querying provenance-backed parameters from YAML."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import yaml

from fermentation.core.tiers import Tier, combine
from fermentation.parameters.schema import Parameter


class ParameterSet:
    """An immutable, name-keyed collection of :class:`Parameter` objects.

    Provides the two things the rest of the engine needs:
      * ``resolve(...)`` — a plain ``{name: float}`` mapping for the hot loop,
      * tier queries — so outputs can declare the confidence of their inputs.
    """

    def __init__(self, parameters: Iterable[Parameter]) -> None:
        store: dict[str, Parameter] = {}
        for p in parameters:
            if p.name in store:
                raise ValueError(f"Duplicate parameter {p.name!r}")
            store[p.name] = p
        self._params = store

    # -- access ---------------------------------------------------------------

    def __getitem__(self, name: str) -> Parameter:
        try:
            return self._params[name]
        except KeyError:
            raise KeyError(f"No parameter {name!r}; have {sorted(self._params)}") from None

    def __contains__(self, name: object) -> bool:
        return name in self._params

    def __iter__(self) -> Any:
        return iter(self._params.values())

    def __len__(self) -> int:
        return len(self._params)

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(self._params)

    def value(self, name: str) -> float:
        return self[name].value

    def tier_of(self, name: str) -> Tier:
        return self[name].tier

    def resolve(self, names: Iterable[str] | None = None) -> dict[str, float]:
        """Return ``{name: value}`` for ``names`` (or all parameters).

        This is what gets handed to :meth:`Process.derivatives`; provenance and
        tier are intentionally dropped here — they belong to setup and analysis,
        not the numeric inner loop.
        """
        if names is None:
            return {n: p.value for n, p in self._params.items()}
        return {n: self[n].value for n in names}

    def lowest_tier(self, names: Iterable[str] | None = None) -> Tier:
        """Lowest (least trustworthy) tier across the named parameters."""
        chosen = self._params if names is None else {n: self[n] for n in names}
        return combine(p.tier for p in chosen.values())

    def tier_map(self) -> dict[str, Tier]:
        """``{name: tier}`` for every parameter — the tier counterpart to
        :meth:`resolve`. Pass it to ``ProcessSet.tier_of``/``tier_map`` (or
        ``simulate(..., param_tiers=...)``) so a Process's output tier is capped by
        the tiers of the parameters it reads (parameter-tier propagation, D-1)."""
        return {n: p.tier for n, p in self._params.items()}

    def merge(self, other: ParameterSet, *, override: bool = False) -> ParameterSet:
        """Combine two sets. By default a name collision is an error; with
        ``override=True``, ``other`` wins (e.g. a strain-specific overlay on top
        of generic defaults)."""
        merged = dict(self._params)
        for name, p in other._params.items():
            if name in merged and not override:
                raise ValueError(f"Parameter {name!r} defined in both sets")
            merged[name] = p
        return ParameterSet(merged.values())


def _parse_mapping(raw: Mapping[str, Any]) -> list[Parameter]:
    params: list[Parameter] = []
    for name, entry in raw.items():
        if not isinstance(entry, Mapping):
            raise ValueError(
                f"Parameter {name!r} must map to a block of fields, got {type(entry).__name__}"
            )
        if "name" in entry and entry["name"] != name:
            raise ValueError(
                f"Parameter key {name!r} disagrees with its 'name' field {entry['name']!r}"
            )
        params.append(Parameter(name=name, **{k: v for k, v in entry.items() if k != "name"}))
    return params


def default_data_dir() -> Path:
    """Filesystem directory holding the packaged parameter YAML files.

    The compile seam locates ``<medium>_<strain>.yaml`` here (and tests read the
    same files). If the package is ever shipped as a wheel that must read its own
    bundled data, switch to ``importlib.resources`` — see the deferred item in
    ``docs/DECISIONS.md``.
    """
    return Path(__file__).resolve().parent / "data"


def load_parameters(*paths: str | Path) -> ParameterSet:
    """Load and validate one or more YAML parameter files into a ``ParameterSet``.

    Each file is a top-level mapping of ``parameter_name -> {value, unit, tier,
    uncertainty, provenance}``. A file (or any entry) missing required provenance
    fields raises ``pydantic.ValidationError`` — load fails loudly rather than
    silently admitting an unsourced magic number. Files are merged left-to-right;
    later files may not redefine an earlier name.
    """
    if not paths:
        raise ValueError("load_parameters() needs at least one path")

    result: ParameterSet | None = None
    for path in paths:
        p = Path(path)
        with p.open("r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
        if raw is None:
            raw = {}
        if not isinstance(raw, Mapping):
            raise ValueError(f"{p}: top level must be a mapping of parameters")
        loaded = ParameterSet(_parse_mapping(raw))
        result = loaded if result is None else result.merge(loaded)
    assert result is not None
    return result
