"""Source registry (data/sources.json) and the adapter registry.

``data/sources.json`` lists every site the project may fetch from, with its
kind, robots.txt snapshot date, terms URL, approval status and rate limit.
An adapter is a Python class registered under the same ``id``. A source runs
only when its registry entry says ``approval_status: "approved"``.
"""

import datetime as _dt
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Type

from horn_drivers.fetch import DEFAULT_CACHE_ROOT, PoliteFetcher

DEFAULT_REGISTRY = Path("data/sources.json")
SOURCE_KINDS = ("manufacturer", "aggregator", "retailer", "reference")
APPROVAL_STATUSES = ("approved", "pending", "excluded", "blocked")
# Floors from the access rules: retailers at least 5 s, loudspeakerdatabase.com 30 s.
MIN_DELAY_FLOOR_S = 5.0


@dataclass
class SourceEntry:
    id: str
    name: str
    kind: str
    base_url: str
    approval_status: str
    min_delay_s: float
    cache_days: float
    robots_url: Optional[str] = None
    robots_snapshot_date: Optional[str] = None
    terms_url: Optional[str] = None
    approval_note: str = ""
    region: Optional[str] = None
    currency: Optional[str] = None
    vat_included: Optional[bool] = None
    status_note: str = ""
    notes: List[str] = field(default_factory=list)

    def check(self) -> List[str]:
        errors = []
        if self.kind not in SOURCE_KINDS:
            errors.append(f"{self.id}: kind must be one of {SOURCE_KINDS}")
        if self.approval_status not in APPROVAL_STATUSES:
            errors.append(f"{self.id}: approval_status must be one of {APPROVAL_STATUSES}")
        if not self.base_url.startswith("https://") and not self.base_url.startswith("http://"):
            errors.append(f"{self.id}: base_url must be http(s)")
        if self.min_delay_s < MIN_DELAY_FLOOR_S:
            errors.append(f"{self.id}: min_delay_s below the {MIN_DELAY_FLOOR_S:g} s floor")
        if self.cache_days <= 0:
            errors.append(f"{self.id}: cache_days must be positive")
        if self.robots_snapshot_date:
            try:
                _dt.date.fromisoformat(self.robots_snapshot_date)
            except ValueError:
                errors.append(f"{self.id}: robots_snapshot_date is not YYYY-MM-DD")
        if self.kind == "retailer" and not (self.region and self.currency):
            errors.append(f"{self.id}: a retailer needs region and currency")
        return errors


def load_source_registry(path: Path = DEFAULT_REGISTRY) -> Dict[str, SourceEntry]:
    data = json.loads(Path(path).read_text())
    entries = {}
    for raw in data["sources"]:
        entry = SourceEntry(**raw)
        problems = entry.check()
        if problems:
            raise ValueError("; ".join(problems))
        if entry.id in entries:
            raise ValueError(f"duplicate source id {entry.id}")
        entries[entry.id] = entry
    return entries


def save_source_registry(entries: Dict[str, SourceEntry], path: Path = DEFAULT_REGISTRY) -> None:
    path = Path(path)
    data = json.loads(path.read_text()) if path.exists() else {}
    data["sources"] = [asdict(entry) for entry in entries.values()]
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def make_fetcher(entry: SourceEntry, cache_root: Path = DEFAULT_CACHE_ROOT, **overrides) -> PoliteFetcher:
    """A PoliteFetcher with the entry's rate limit and cache window."""
    return PoliteFetcher(source_id=entry.id, min_delay_s=entry.min_delay_s,
                         cache_days=entry.cache_days, cache_root=Path(cache_root), **overrides)


# -- adapter classes -------------------------------------------------------------

ADAPTERS: Dict[str, Type] = {}


def register_adapter(cls):
    """Class decorator: register a SourceAdapter subclass under its source_id."""
    source_id = getattr(cls, "source_id", None)
    if not source_id:
        raise ValueError(f"{cls.__name__} has no source_id")
    if source_id in ADAPTERS and ADAPTERS[source_id] is not cls:
        raise ValueError(f"adapter for {source_id} already registered")
    ADAPTERS[source_id] = cls
    return cls


def adapter_ids() -> List[str]:
    _import_builtin_adapters()
    return sorted(ADAPTERS)


def get_adapter_class(source_id: str):
    _import_builtin_adapters()
    if source_id not in ADAPTERS:
        raise KeyError(f"no adapter registered for {source_id}; known: {', '.join(sorted(ADAPTERS)) or 'none'}")
    return ADAPTERS[source_id]


def create_adapter(source_id: str, registry_path: Path = DEFAULT_REGISTRY,
                   cache_root: Path = DEFAULT_CACHE_ROOT, fetcher: Optional[PoliteFetcher] = None):
    """Adapter instance for an approved source, wired to a PoliteFetcher."""
    entries = load_source_registry(registry_path)
    if source_id not in entries:
        raise KeyError(f"{source_id} is not in {registry_path}")
    entry = entries[source_id]
    if entry.approval_status != "approved":
        raise PermissionError(f"{source_id} has approval_status {entry.approval_status!r}: {entry.approval_note}")
    cls = get_adapter_class(source_id)
    if cls.kind != entry.kind:
        raise ValueError(f"adapter kind {cls.kind} does not match registry kind {entry.kind}")
    return cls(fetcher or make_fetcher(entry, cache_root), entry)


def _import_builtin_adapters() -> None:
    """Import every module in horn_drivers.sources so their decorators run."""
    import importlib
    import pkgutil

    import horn_drivers.sources as package
    for module in pkgutil.iter_modules(package.__path__):
        if module.name not in ("base", "registry", "cli"):
            importlib.import_module(f"horn_drivers.sources.{module.name}")
