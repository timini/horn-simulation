"""Driver and price source adapters.

See ``base`` for the interface and ``registry`` for data/sources.json and
adapter registration. Each adapter lives in its own module in this package
and registers itself with ``@register_adapter``.
"""

from horn_drivers.sources.base import (DriverIndex, DriverResult, DriverSource, Listing, PriceSource,
                                       ProductRef, RunSummary, SourceAdapter)
from horn_drivers.sources.registry import (ADAPTERS, SourceEntry, adapter_ids, create_adapter, get_adapter_class,
                                           load_source_registry, make_fetcher, register_adapter,
                                           save_source_registry)

__all__ = ["DriverIndex", "DriverResult", "DriverSource", "Listing", "PriceSource", "ProductRef", "RunSummary",
           "SourceAdapter", "ADAPTERS", "SourceEntry", "adapter_ids", "create_adapter", "get_adapter_class",
           "load_source_registry", "make_fetcher", "register_adapter", "save_source_registry"]
