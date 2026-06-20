"""The state vector and its schema.

Design decision (see ``docs/DECISIONS.md``): the *integrated* state is a plain,
contiguous ``float64`` numpy array so that ``scipy.integrate.solve_ivp`` can drive
it with no per-element Python objects in the hot loop. The mapping from physical
variable names to array indices — including which variables are vectors (e.g. the
multi-sugar ``S`` for beer) — lives here in :class:`StateSchema`. Tier and
provenance metadata deliberately do **not** ride inside these floats; they are
derived at the analysis boundary from the Processes and parameters that fed each
variable.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class VarSpec:
    """Specification of one state variable.

    A variable occupies ``size`` contiguous slots in the flat state array.
    ``size == 1`` is a scalar (biomass, ethanol, temperature, …); ``size > 1``
    is a vector, used for sequentially-consumed beer sugars (glucose, maltose,
    maltotriose) where ``components`` names each slot.
    """

    name: str
    unit: str
    size: int = 1
    description: str = ""
    components: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.size < 1:
            raise ValueError(f"VarSpec {self.name!r}: size must be >= 1, got {self.size}")
        if self.components and len(self.components) != self.size:
            raise ValueError(
                f"VarSpec {self.name!r}: {len(self.components)} component names "
                f"for size {self.size}"
            )


class StateSchema:
    """An ordered registry mapping variable names to slices of the state array.

    The schema is immutable once constructed. It is the single source of truth
    for how a physical state maps onto the flat numpy vector the solver sees.
    """

    def __init__(self, specs: Sequence[VarSpec]) -> None:
        if not specs:
            raise ValueError("StateSchema requires at least one variable")
        names = [s.name for s in specs]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate variable names in schema: {names}")

        self._specs: tuple[VarSpec, ...] = tuple(specs)
        self._slices: dict[str, slice] = {}
        cursor = 0
        for spec in self._specs:
            self._slices[spec.name] = slice(cursor, cursor + spec.size)
            cursor += spec.size
        self._size = cursor

    @property
    def size(self) -> int:
        """Total length of the flat state array."""
        return self._size

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(s.name for s in self._specs)

    @property
    def specs(self) -> tuple[VarSpec, ...]:
        return self._specs

    def spec(self, name: str) -> VarSpec:
        for s in self._specs:
            if s.name == name:
                return s
        raise KeyError(f"Unknown state variable {name!r}; known: {self.names}")

    def slice(self, name: str) -> slice:
        """Slice into the flat array occupied by ``name``."""
        try:
            return self._slices[name]
        except KeyError:
            raise KeyError(f"Unknown state variable {name!r}; known: {self.names}") from None

    def pack(self, values: Mapping[str, float | Sequence[float]]) -> FloatArray:
        """Build a flat state array from a name → value(s) mapping.

        Every variable must be supplied exactly once; scalar variables take a
        number, vector variables take a sequence of the right length.
        """
        missing = set(self.names) - set(values)
        extra = set(values) - set(self.names)
        if missing:
            raise ValueError(f"pack() missing values for: {sorted(missing)}")
        if extra:
            raise ValueError(f"pack() got unknown variables: {sorted(extra)}")

        arr = np.empty(self._size, dtype=np.float64)
        for spec in self._specs:
            block = np.atleast_1d(np.asarray(values[spec.name], dtype=np.float64))
            if block.size != spec.size:
                raise ValueError(
                    f"Variable {spec.name!r} expects {spec.size} value(s), got {block.size}"
                )
            arr[self._slices[spec.name]] = block
        return arr

    def unpack(self, array: FloatArray) -> dict[str, float | FloatArray]:
        """Inverse of :meth:`pack` — scalar vars come back as floats."""
        self._check_shape(array)
        out: dict[str, float | FloatArray] = {}
        for spec in self._specs:
            block = array[self._slices[spec.name]]
            out[spec.name] = float(block[0]) if spec.size == 1 else block.copy()
        return out

    def get(self, array: FloatArray, name: str) -> float | FloatArray:
        """Read one variable out of a flat array."""
        self._check_shape(array)
        spec = self.spec(name)
        block = array[self._slices[name]]
        return float(block[0]) if spec.size == 1 else block.copy()

    def zeros(self) -> FloatArray:
        """A zero-filled state array matching this schema."""
        return np.zeros(self._size, dtype=np.float64)

    def _check_shape(self, array: FloatArray) -> None:
        if array.shape != (self._size,):
            raise ValueError(f"State array has shape {array.shape}, expected ({self._size},)")

    def __len__(self) -> int:
        return self._size

    def __contains__(self, name: object) -> bool:
        return name in self._slices

    def __repr__(self) -> str:
        parts = ", ".join(f"{s.name}[{s.size}]" if s.size > 1 else s.name for s in self._specs)
        return f"StateSchema({parts}; size={self._size})"


@dataclass
class StateVector:
    """Ergonomic schema + array pairing for use *outside* the integration loop.

    The solver works on raw arrays; this convenience view is for setting up
    initial conditions, inspecting results, and writing readable tests.
    """

    schema: StateSchema
    array: FloatArray = field(default_factory=lambda: np.empty(0))

    def __post_init__(self) -> None:
        if self.array.size == 0:
            self.array = self.schema.zeros()
        else:
            self.schema._check_shape(self.array)

    @classmethod
    def from_values(
        cls, schema: StateSchema, values: Mapping[str, float | Sequence[float]]
    ) -> StateVector:
        return cls(schema=schema, array=schema.pack(values))

    def __getitem__(self, name: str) -> float | FloatArray:
        return self.schema.get(self.array, name)

    def __setitem__(self, name: str, value: float | Sequence[float]) -> None:
        block = np.atleast_1d(np.asarray(value, dtype=np.float64))
        sl = self.schema.slice(name)
        if block.size != (sl.stop - sl.start):
            raise ValueError(f"Variable {name!r} expects {sl.stop - sl.start} value(s)")
        self.array[sl] = block

    def as_dict(self) -> dict[str, float | FloatArray]:
        return self.schema.unpack(self.array)
