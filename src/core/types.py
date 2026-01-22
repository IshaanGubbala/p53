from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence, Union

Mutation = str
MutationSet = tuple[Mutation, ...]
JsonLike = dict[str, Any]
PathLike = Union[str, Path]

# Generic containers used across the pipeline.
Record = Mapping[str, Any]
Records = Sequence[Record]
