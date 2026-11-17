"""driver-base: an open driver dataset (CC BY 4.0) compiled from manufacturer pages.

driver-base (https://github.com/JonathanDCohen/driver-base) scrapes the
product pages of 15 manufacturers and publishes one JSON file,
``web/drivers.json``. The dataset is licensed CC BY 4.0, which allows reuse
with attribution and a note of what was changed (see ``ATTRIBUTION``). The
import is pinned to one commit so a rerun gives the same records.

What the import does:

- keeps only values that driver-base read from a page or PDF. Values it
  derived or inferred (for example program power as twice AES power, or
  1 W sensitivity computed from 2.83 V) are left out;
- converts units to SI and builds version 2 records, with a field_provenance
  entry per field that names the manufacturer page driver-base read, the
  date it read it, and the dataset it came through;
- matches each record against the existing database by canonical
  manufacturer, normalised model and impedance. A matched record is never
  overwritten. Values that disagree are written to the conflict report. A
  matched scraped record (``data/drivers``, version 1) that lacks power,
  sensitivity or band gets those fields added, but only when Fs, Re, Bl and
  Sd agree, so both records clearly describe the same driver;
- checks each record for internal consistency: Fs against Mms and Cms
  (10 %), Qes against Fs, Mms, Re and Bl (20 %), and the page's Vas
  against Sd, Fs and Mms (a factor of 1.5). A Cms that disagrees
  while Qes agrees is dropped; any other disagreement sends the record to
  the incomplete folder with the reason;
- copies driver-base's warn_flags into the record notes and marks a power
  rating it flags as not AES in that field's provenance;
- writes compression drivers without a full Thiele/Small set to the
  incomplete folder, and skips passive radiators, shakers and makers not on
  the canonical list.
"""

import csv
import json
import math
import re
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from horn_drivers.fetch import FetchResult
from horn_drivers.identity import canonical_manufacturer, match_key, record_identity
from horn_drivers.records import check_structure, driver_id_for, make_record, provenance_entry, write_incomplete, write_record
from horn_drivers.sources.base import DriverResult, DriverSource, ProductRef, RunSummary
from horn_drivers.sources.registry import register_adapter

DATASET_REPO = "JonathanDCohen/driver-base"
DATASET_COMMIT = "f0fc602b2fb4154b20b90d4cae03fa82296b249a"
DATASET_PATH = "web/drivers.json"
DATASET_URL = f"https://raw.githubusercontent.com/{DATASET_REPO}/{DATASET_COMMIT}/{DATASET_PATH}"
LICENSE = "CC BY 4.0"
ATTRIBUTION = ("Data from driver-base (https://github.com/JonathanDCohen/driver-base), "
               "licensed under CC BY 4.0. Changed: units converted to SI, derived values dropped, "
               "records matched to this database.")

# driver-base spec_source values that mean the number was read from the source.
PAGE_SOURCES = frozenset({"html_table", "html_div_pairs", "html_prose", "html_grid", "inline_js", "json_api"})
PDF_SOURCES = frozenset({"pdf_table", "pdf_text"})

# driver-base field, our parameter, scale to SI
FIELD_MAP: Tuple[Tuple[str, str, float], ...] = (
    ("fs_hz", "fs_hz", 1.0),
    ("re_ohm", "re_ohm", 1.0),
    ("bl_tm", "bl_tm", 1.0),
    ("sd_cm2", "sd_m2", 1e-4),
    ("mms_g", "mms_kg", 1e-3),
    ("le_mh", "le_h", 1e-3),
    ("qms", "qms", 1.0),
    ("qes", "qes", 1.0),
    ("qts", "qts", 1.0),
    ("xmax_mm", "xmax_m", 1e-3),
    ("cms_mm_per_n", "cms_m_per_n", 1e-3),
    ("rms_ns_per_m", "rms_kg_per_s", 1.0),
    ("impedance_nominal_ohm", "nominal_impedance_ohm", 1.0),
    ("power_aes_watts", "power_w", 1.0),
    ("power_program_watts", "peak_power_w", 1.0),
    ("sensitivity_db_1w_1m", "sensitivity_db_1w1m", 1.0),
)
BAND_MAP = (("freq_low_hz", "usable_f_low_hz"), ("freq_high_hz", "usable_f_high_hz"))
ESSENTIAL = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg", "le_h")

DRIVER_TYPES = {"lf_woofer": "cone", "guitar_bass": "cone", "fullrange": "cone", "coax": "cone",
                "hf_compression": "compression"}
SKIPPED_KINDS = {"passive": "passive radiator (no motor)", "shaker": "tactile transducer",
                 "tweeter": "dome tweeter (not a horn driver here)"}

# Fields compared between a driver-base record and an existing one, with the
# relative tolerance from the plan's conflict rule (5 %, 10 % for Le and Xmax).
COMPARE_TOLERANCE = {"fs_hz": 0.05, "re_ohm": 0.05, "bl_tm": 0.05, "sd_m2": 0.05, "mms_kg": 0.05,
                     "le_h": 0.10, "qts": 0.05, "qes": 0.05, "qms": 0.05, "xmax_m": 0.10,
                     "power_w": 0.05, "sensitivity_db_1w1m": 0.01}
IDENTITY_FIELDS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2")
# A record whose model name differs only by a suffix (an impedance or coil
# code such as "HSB322-22") and whose Fs, Re, Bl, Sd and Mms all agree within
# this tolerance is the same driver.
SAME_VALUES_TOLERANCE = 0.01
SAME_VALUES_FIELDS = ("fs_hz", "re_ohm", "bl_tm", "sd_m2", "mms_kg")
FILL_FIELDS = ("power_w", "peak_power_w", "sensitivity_db_1w1m")
# Internal consistency gate: Fs = 1 / (2 pi sqrt(Mms Cms)) and
# Qes = 2 pi Fs Mms Re / Bl^2, each within this relative tolerance.
FS_CONSISTENCY_TOLERANCE = 0.10
QES_CONSISTENCY_TOLERANCE = 0.20
# Vas = rho c^2 Sd^2 Cms, with Cms from Fs and Mms: the page's Vas (used only
# for this check, not stored) must agree within this factor. It catches Sd or
# Mms parse slips by a factor of ten.
VAS_CONSISTENCY_FACTOR = 1.5
RHO_C2 = 1.204 * 343.0 ** 2
# A page nominal size is used when the effective cone diameter from Sd is
# this fraction of it; otherwise the size is guessed from Sd.
NOMINAL_SIZE_SD_RATIO = (0.6, 1.0)
NON_AES_POWER_FLAG = "power_std_non_aes"
# Nominal sizes (inches) used when driver-base read the size from the page.
NOMINAL_SIZES_IN = (2, 3, 4, 5, 6, 6.5, 7, 8, 10, 12, 15, 18, 21, 24)
PREFIX_REVIEW_COLUMNS = ("canonical_id", "driver_base_model", "existing_driver_id", "path", "extra_suffix",
                         "driver_base_url", "reason")
CONFLICT_COLUMNS = ("driver_id", "path", "field", "existing", "driver_base", "relative_difference",
                    "tolerance", "driver_base_url", "note")


def _number(value) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) and value > 0 else None


def _method(spec_source: Optional[str]) -> Optional[str]:
    if spec_source in PAGE_SOURCES:
        return "page"
    if spec_source in PDF_SOURCES:
        return "pdf"
    return None  # derived, inferred, override or missing: not imported


def consistency_problems(parameters: Dict[str, float], vas_m3: Optional[float] = None) -> Dict[str, str]:
    """Checks that fail, as {"fs_cms" | "qes" | "vas": description}. Checks whose inputs are missing are skipped."""
    problems = {}
    fs, mms, cms = parameters.get("fs_hz"), parameters.get("mms_kg"), parameters.get("cms_m_per_n")
    if fs and mms and cms:
        computed = 1.0 / (2.0 * math.pi * math.sqrt(mms * cms))
        if abs(computed - fs) / fs > FS_CONSISTENCY_TOLERANCE:
            problems["fs_cms"] = f"Fs {fs:g} Hz but Mms and Cms give {computed:.3g} Hz"
    re_ohm, bl, qes = parameters.get("re_ohm"), parameters.get("bl_tm"), parameters.get("qes")
    if fs and mms and re_ohm and bl and qes:
        computed = 2.0 * math.pi * fs * mms * re_ohm / bl ** 2
        if abs(computed - qes) / qes > QES_CONSISTENCY_TOLERANCE:
            problems["qes"] = f"Qes {qes:g} but Fs, Mms, Re and Bl give {computed:.3g}"
    sd = parameters.get("sd_m2")
    if fs and mms and sd and vas_m3:
        computed = RHO_C2 * sd ** 2 / ((2.0 * math.pi * fs) ** 2 * mms)
        if max(computed / vas_m3, vas_m3 / computed) > VAS_CONSISTENCY_FACTOR:
            problems["vas"] = f"page Vas {vas_m3 * 1e3:.3g} L but Sd, Fs and Mms give {computed * 1e3:.3g} L"
    return problems


def nominal_size_label(size_mm: Optional[float], sd_m2: Optional[float] = None) -> Optional[str]:
    """The nearest standard nominal size, as "8in" or "6.5in", for a size in mm.

    None when the size is missing, or when Sd is given and its effective
    diameter is not a plausible fraction of the size (a misread size).
    """
    if size_mm is None:
        return None
    if sd_m2 is not None:
        ratio = 2.0 * math.sqrt(sd_m2 / math.pi) * 1e3 / size_mm
        if not NOMINAL_SIZE_SD_RATIO[0] <= ratio <= NOMINAL_SIZE_SD_RATIO[1]:
            return None
    inches = size_mm / 25.4
    best = min(NOMINAL_SIZES_IN, key=lambda size: abs(size - inches))
    return f"{best:g}in"


def convert_entry(entry: dict, *, impedance_in_id: bool = False, dataset_url: str = DATASET_URL) -> DriverResult:
    """One driver-base record as a version 2 record, or a reason it is skipped.

    Pure: no network, no file access. ``impedance_in_id`` appends the
    impedance to driver_id, for models listed in several impedances.
    """
    ref = ProductRef(url=(entry.get("source_urls") or [dataset_url])[0], manufacturer=entry.get("manufacturer"),
                     model=entry.get("model"), hints=(("canonical_id", str(entry.get("canonical_id"))),))
    kind = entry.get("driver_kind")
    if kind in SKIPPED_KINDS:
        return DriverResult(ref, None, f"skipped kind: {SKIPPED_KINDS[kind]}")
    manufacturer = canonical_manufacturer(entry.get("manufacturer"))
    if manufacturer is None:
        return DriverResult(ref, None, f"manufacturer {entry.get('manufacturer')!r} is not on the canonical list")
    model = (entry.get("model") or "").strip()
    if not model:
        return DriverResult(ref, None, "no model name")
    urls = [u for u in entry.get("source_urls") or [] if str(u).startswith(("http://", "https://"))]
    if not urls:
        return DriverResult(ref, None, "no source URL")
    url = urls[0]
    read_on = str(entry.get("fetched_at") or entry.get("last_scraped_at") or "")[:10]
    sources = entry.get("spec_source") or {}

    parameters: Dict[str, float] = {}
    provenance: Dict[str, dict] = {}
    dropped: List[str] = []

    def take(src: str, dest: str, scale: float):
        value = _number(entry.get(src))
        if value is None:
            return
        method = _method(sources.get(src))
        if method is None:
            dropped.append(src)
            return
        parameters[dest] = float(f"{value * scale:.12g}")
        provenance[dest] = provenance_entry("driver-base", url, read_on, method, via=dataset_url,
                                            upstream_field=src, upstream_spec_source=sources.get(src))

    for src, dest, scale in FIELD_MAP:
        take(src, dest, scale)
    band = {}
    for src, dest in BAND_MAP:
        value = _number(entry.get(src))
        method = _method(sources.get(src))
        if value is not None and method is not None:
            band[dest] = value
            provenance[dest] = provenance_entry("driver-base", url, read_on, method, via=dataset_url,
                                                upstream_field=src, upstream_spec_source=sources.get(src))
    if len(band) != 2 or band["usable_f_low_hz"] >= band["usable_f_high_hz"]:
        for dest in band:
            provenance.pop(dest, None)
        band = {}
    # A peak rating below the nominal one is a parse slip, not a rating.
    if "peak_power_w" in parameters and parameters["peak_power_w"] < parameters.get("power_w", 0.0):
        parameters.pop("peak_power_w")
        provenance.pop("peak_power_w")
    warn_flags = [str(flag) for flag in entry.get("warn_flags") or []]
    if "power_w" in provenance and NON_AES_POWER_FLAG in warn_flags:
        provenance["power_w"]["power_standard"] = "not AES (driver-base warn flag power_std_non_aes)"

    vas_l = _number(entry.get("vas_liters")) if _method(sources.get("vas_liters")) else None
    problems = consistency_problems(parameters, vas_l * 1e-3 if vas_l else None)
    inconsistent = None
    dropped_cms = None
    if "fs_cms" in problems and "qes" not in problems and "qes" in parameters:
        # Mms agrees with Fs, Re, Bl and Qes, so the stated Cms is the bad value.
        dropped_cms = problems.pop("fs_cms")
        parameters.pop("cms_m_per_n")
        provenance.pop("cms_m_per_n")
    if problems:
        inconsistent = "inconsistent Thiele/Small values: " + "; ".join(problems[k] for k in sorted(problems))

    driver_type = DRIVER_TYPES.get(kind, "unknown")
    type_source = (f"driver-base driver_kind {kind!r}, from the manufacturer's product category"
                   if driver_type != "unknown" else None)
    impedance = parameters.get("nominal_impedance_ohm")
    notes = [f"Imported from driver-base commit {DATASET_COMMIT[:12]} ({LICENSE}); "
             f"driver-base kind {kind!r}, canonical_id {entry.get('canonical_id')!r}."]
    if kind == "coax":
        notes.append("Coaxial driver: the parameters are for the low-frequency section.")
    if dropped:
        notes.append("Derived or inferred driver-base values not imported: " + ", ".join(sorted(dropped)) + ".")
    if dropped_cms:
        notes.append(f"Cms not imported: {dropped_cms}, while Qes agrees with Mms.")
    if warn_flags:
        notes.append("driver-base warn_flags: " + ", ".join(warn_flags) + ".")
    nominal_diameter = None
    if driver_type == "cone":
        if _method((entry.get("spec_source") or {}).get("nominal_size_mm")) is not None:
            nominal_diameter = nominal_size_label(_number(entry.get("nominal_size_mm")), parameters.get("sd_m2"))
        if nominal_diameter is None and "sd_m2" in parameters:
            from horn_drivers.scraper import infer_nominal_diameter
            nominal_diameter = infer_nominal_diameter(parameters["sd_m2"])
    try:
        record = make_record(manufacturer=manufacturer, model=model, parameters=parameters,
                             field_provenance=provenance, driver_type=driver_type, impedance_ohm=impedance,
                             nominal_diameter=nominal_diameter, usable_f_low_hz=band.get("usable_f_low_hz"),
                             usable_f_high_hz=band.get("usable_f_high_hz"), id_includes_impedance=impedance_in_id,
                             notes=notes)
    except ValueError as error:
        return DriverResult(ref, None, f"record not built: {error}")
    if type_source:
        record["type_source"] = type_source
    missing = [name for name in ESSENTIAL if name not in parameters]
    if not missing and "qts" not in parameters and not ("qes" in parameters and "qms" in parameters):
        missing = ["qts"]
    if missing:
        if driver_type == "compression":
            return DriverResult(ref, record, "missing Thiele/Small fields: " + ", ".join(missing), incomplete=True)
        return DriverResult(ref, None, "missing Thiele/Small fields: " + ", ".join(missing))
    if inconsistent:
        return DriverResult(ref, record, inconsistent, incomplete=True)
    return DriverResult(ref, record)


def impedance_variants(entries: Sequence[dict]) -> Dict[Tuple[str, str], int]:
    """How many impedance variants each (manufacturer, model) has in the dataset."""
    seen: Dict[Tuple[str, str], set] = {}
    for entry in entries:
        manufacturer = canonical_manufacturer(entry.get("manufacturer"))
        if manufacturer and entry.get("model"):
            key = match_key(manufacturer, entry["model"], None)[:2]
            seen.setdefault(key, set()).add(entry.get("impedance_nominal_ohm"))
    return {key: len(values) for key, values in seen.items()}


def shared_keys(entries: Sequence[dict]) -> Dict[Tuple[str, str, Optional[float]], int]:
    """Match keys that more than one dataset record maps to, with their counts.

    driver-base sometimes gives two products the same model string (for
    example two Monacor SP-6 versions); neither can be matched safely.
    """
    counts: Dict[Tuple, int] = {}
    for entry in entries:
        manufacturer = canonical_manufacturer(entry.get("manufacturer"))
        if manufacturer and entry.get("model") and entry.get("driver_kind") not in SKIPPED_KINDS:
            key = match_key(manufacturer, entry["model"], _number(entry.get("impedance_nominal_ohm")))
            counts[key] = counts.get(key, 0) + 1
    return {key: count for key, count in counts.items() if count > 1}


def record_key(record: dict) -> Optional[Tuple[str, str, Optional[float]]]:
    identity = record_identity(record)
    try:
        return match_key(identity.get("manufacturer") or "", identity.get("model") or "",
                         (identity.get("variant") or {}).get("impedance_ohm")
                         or (record.get("parameters") or {}).get("nominal_impedance_ohm"))
    except ValueError:
        return None


def model_suffix(existing: dict, new: dict) -> Optional[str]:
    """The text one normalised model name adds to the other, or None when neither extends the other."""
    old_key, new_key = record_key(existing), record_key(new)
    if old_key is None or new_key is None or old_key[0] != new_key[0]:
        return None
    longer, shorter = sorted((old_key[1], new_key[1]), key=len, reverse=True)
    return longer[len(shorter):] if longer.startswith(shorter) else None


def is_impedance_suffix(suffix: str, impedances: Sequence[Optional[float]]) -> bool:
    """True for "" or an impedance or dual-coil code ("4", "22", "8OHM") naming one of *impedances*."""
    if suffix == "":
        return True
    match = re.fullmatch(r"(\d{1,2})(OHM)?", suffix)
    if not match:
        return False
    digits = match.group(1)
    for impedance in impedances:
        if impedance and float(impedance).is_integer():
            code = str(int(impedance))
            if digits in (code, code * 2):
                return True
    return False


def _impedances(record: dict) -> List[Optional[float]]:
    variant = (record_identity(record).get("variant") or {})
    return [variant.get("impedance_ohm"), (record.get("parameters") or {}).get("nominal_impedance_ohm")]


def same_values(existing: dict, new: dict) -> bool:
    """Fs, Re, Bl, Sd and Mms all present and within SAME_VALUES_TOLERANCE."""
    for name in SAME_VALUES_FIELDS:
        old, value = _value(existing, name), _value(new, name)
        if old is None or value is None or abs(value - old) / old > SAME_VALUES_TOLERANCE:
            return False
    return True


def same_driver_by_values(existing: dict, new: dict) -> bool:
    """True when the model names differ only by an impedance or coil code and the core values agree.

    Catches the same driver filed under a coil or impedance suffix ("HSB322"
    and "HSB322-22") or under a different nominal impedance convention for a
    dual-coil driver. A suffix that is not an impedance code ("EVO", "WP",
    "-2" on an 8 ohm driver) may name a different revision, so it does not
    match; see prefix_candidate.
    """
    suffix = model_suffix(existing, new)
    if suffix is None or not is_impedance_suffix(suffix, _impedances(existing) + _impedances(new)):
        return False
    return same_values(existing, new)


def prefix_candidate(existing: dict, new: dict) -> Optional[str]:
    """The suffix when the names differ by a non-impedance suffix but the core values agree, else None."""
    suffix = model_suffix(existing, new)
    if suffix is None or is_impedance_suffix(suffix, _impedances(existing) + _impedances(new)):
        return None
    return suffix if same_values(existing, new) else None


def _value(record: dict, name: str) -> Optional[float]:
    if name in ("usable_f_low_hz", "usable_f_high_hz"):
        return _number(record.get(name))
    return _number((record.get("parameters") or {}).get(name))


def compare_records(existing: dict, new: dict) -> List[dict]:
    """Fields where the two records disagree beyond the conflict tolerance."""
    rows = []
    for name, tolerance in COMPARE_TOLERANCE.items():
        old, value = _value(existing, name), _value(new, name)
        if old is None or value is None:
            continue
        difference = abs(value - old) / old
        if difference > tolerance:
            rows.append({"field": name, "existing": old, "driver_base": value,
                         "relative_difference": round(difference, 4), "tolerance": tolerance})
    return rows


def fill_missing(existing: dict, new: dict) -> Tuple[Optional[dict], List[str]]:
    """Copy power, sensitivity and band from *new* into a scraped record that lacks them.

    Returns (updated record, filled field names), or (None, []) when nothing
    is filled. Fills only when Fs, Re, Bl and Sd are all present on both
    sides and agree within tolerance. Existing values are never changed.
    """
    for name in IDENTITY_FIELDS:
        old, value = _value(existing, name), _value(new, name)
        if old is None or value is None or abs(value - old) / old > COMPARE_TOLERANCE[name]:
            return None, []
    params = dict(existing.get("parameters") or {})
    provenance = dict(existing.get("field_provenance") or {})
    new_provenance = new.get("field_provenance") or {}
    filled = []
    for name in FILL_FIELDS:
        if params.get(name) is None and _value(new, name) is not None:
            params[name] = new["parameters"][name]
            provenance[name] = new_provenance[name]
            filled.append(name)
    # A filled peak rating below the filled nominal one is dropped; an
    # existing peak rating is never touched.
    if "power_w" in filled and "peak_power_w" in filled and params["peak_power_w"] < params["power_w"]:
        params.pop("peak_power_w")
        provenance.pop("peak_power_w")
        filled.remove("peak_power_w")
    updated = {**existing, "parameters": params}
    if existing.get("usable_f_low_hz") is None and existing.get("usable_f_high_hz") is None \
            and new.get("usable_f_low_hz") is not None:
        updated["usable_f_low_hz"], updated["usable_f_high_hz"] = new["usable_f_low_hz"], new["usable_f_high_hz"]
        provenance["usable_f_low_hz"] = new_provenance["usable_f_low_hz"]
        provenance["usable_f_high_hz"] = new_provenance["usable_f_high_hz"]
        filled += ["usable_f_low_hz", "usable_f_high_hz"]
    if not filled:
        return None, []
    updated["field_provenance"] = provenance
    meta = dict(updated.get("provenance") or {})
    limitations = [note for note in meta.get("known_limitations", [])
                   if not note.startswith("Band, power handling, Rms and sensitivity were not captured")]
    if "known_limitations" in meta:
        still_missing = [label for label, name in (("band", "usable_f_low_hz"), ("power handling", "power_w"),
                                                   ("sensitivity", "sensitivity_db_1w1m"))
                         if name not in filled and _value(updated, name) is None] + ["Rms"]
        meta["known_limitations"] = limitations + [", ".join(still_missing).capitalize() + " not captured."]
    meta["supplemented_fields"] = sorted(set(meta.get("supplemented_fields", [])) | set(filled))
    meta["supplemented_from"] = ("driver-base (CC BY 4.0) values read from the manufacturer page; "
                                 "see field_provenance for each field.")
    updated["provenance"] = meta
    return updated, filled


class ImportSummary(RunSummary):
    def __init__(self, source_id: str):
        super().__init__(source_id)
        self.matched = 0
        self.filled = 0
        self.conflicts = 0
        self.duplicates = 0
        self.prefix_review = 0
        self.skip_reasons: Dict[str, int] = {}

    def to_dict(self) -> dict:
        return {**super().to_dict(), "matched_existing": self.matched, "filled_existing": self.filled,
                "conflicting_fields": self.conflicts, "ambiguous_in_dataset": self.duplicates,
                "prefix_matches_for_review": self.prefix_review,
                "skip_reasons": dict(sorted(self.skip_reasons.items(), key=lambda item: -item[1]))}


def _skip_reason_class(reason: str) -> str:
    return reason.split(":")[0] if reason.startswith(("missing", "skipped kind", "record not built")) else reason


@register_adapter
class DriverBaseSource(DriverSource):
    source_id = "driver-base"
    kind = "aggregator"

    def __init__(self, fetcher, entry=None):
        super().__init__(fetcher, entry)
        self._parsed: Dict[str, dict] = {}

    def dataset(self, page: FetchResult) -> dict:
        """The parsed dataset, parsed once per run."""
        key = page.url
        if key not in self._parsed:
            data = json.loads(page.text)
            self._parsed[key] = {"meta": {k: v for k, v in data.items() if k != "drivers"},
                                 "drivers": data["drivers"],
                                 "by_id": {d["canonical_id"]: d for d in reversed(data["drivers"])},
                                 "variants": impedance_variants(data["drivers"]),
                                 "shared": shared_keys(data["drivers"])}
        return self._parsed[key]

    def discover(self, limit: Optional[int] = None) -> Iterator[ProductRef]:
        page = self.fetcher.get(DATASET_URL)
        if not page.ok:
            return
        for count, entry in enumerate(self.dataset(page)["drivers"]):
            if limit is not None and count >= limit:
                return
            yield ProductRef(url=DATASET_URL, manufacturer=entry.get("manufacturer"), model=entry.get("model"),
                             hints=(("canonical_id", entry["canonical_id"]), ("index", str(count))))

    def parse_driver(self, page: FetchResult, ref: ProductRef) -> DriverResult:
        data = self.dataset(page)
        entry = data["drivers"][int(ref.hint("index"))]
        manufacturer = canonical_manufacturer(entry.get("manufacturer"))
        several = bool(manufacturer and entry.get("model")
                       and data["variants"].get(match_key(manufacturer, entry["model"], None)[:2], 0) > 1)
        return convert_entry(entry, impedance_in_id=several, dataset_url=page.url)

    def run(self, db_dir: Path, incomplete_dir: Path, limit: Optional[int] = None,
            existing_dirs: Optional[Sequence[Path]] = None, report_dir: Optional[Path] = None,
            fill: bool = True) -> ImportSummary:
        """Import new records, fill matched scraped records, report conflicts.

        *existing_dirs* (default: *db_dir* and ``data/drivers-curated``) are
        the folders deduplicated against. Only records under *db_dir* with
        record_schema_version 1 are ever filled. Reports go to *report_dir*
        (default ``data/sources/driver-base``).
        """
        from horn_drivers.loader import _driver_from_dict
        db_dir, incomplete_dir = Path(db_dir), Path(incomplete_dir)
        report_dir = Path(report_dir) if report_dir else Path("data/sources/driver-base")
        dirs = [Path(d) for d in (existing_dirs or [db_dir, Path("data/drivers-curated")])]
        existing: Dict[Tuple, Tuple[Path, dict]] = {}
        by_maker: Dict[str, List[Tuple[Path, dict]]] = {}
        taken_paths = set()
        for folder in dirs:
            for path in sorted(folder.glob("*/*.json")) if folder.is_dir() else []:
                try:
                    record = json.loads(path.read_text())
                except (json.JSONDecodeError, OSError):
                    continue
                taken_paths.add(path.resolve())
                key = record_key(record)
                if key is not None:
                    existing.setdefault(key, (path, record))
                    by_maker.setdefault(key[0], []).append((path, record))

        summary = ImportSummary(self.source_id)
        conflicts: List[dict] = []
        matches: List[dict] = []
        review: List[dict] = []
        imported: Dict[Tuple, str] = {}
        for ref, page in self._pages(limit, summary):
            result = self.parse_driver(page, ref)
            shared = self.dataset(page)["shared"]
            if result.record is None:
                reason = result.reason or "no record"
                summary.skipped.append((ref.hint("canonical_id"), reason))
                cls = _skip_reason_class(reason)
                summary.skip_reasons[cls] = summary.skip_reasons.get(cls, 0) + 1
                continue
            record = result.record
            key = record_key(record)
            if key in shared:
                summary.duplicates += 1
                summary.skipped.append((ref.hint("canonical_id"),
                                        f"ambiguous: {shared[key]} dataset records share this manufacturer, model and impedance"))
                continue
            hit = existing.get(key) or (existing.get(key[:2] + (None,)) if key[2] is not None else None)
            how = "manufacturer, model and impedance"
            if hit is None:
                hit = next((item for item in by_maker.get(key[0], []) if same_driver_by_values(item[1], record)), None)
                how = "impedance or coil suffix and Fs, Re, Bl, Sd, Mms within 1 %"
            if hit is None:
                candidate = next(((item, suffix) for item in by_maker.get(key[0], [])
                                  if (suffix := prefix_candidate(item[1], record)) is not None), None)
                if candidate is not None:
                    (path, old), suffix = candidate
                    review.append({"canonical_id": ref.hint("canonical_id"), "driver_base_model": record["model_name"],
                                   "existing_driver_id": old.get("driver_id"), "path": str(path),
                                   "extra_suffix": suffix, "driver_base_url": ref.url,
                                   "reason": "model names differ by a suffix that is not an impedance code; "
                                             "Fs, Re, Bl, Sd and Mms agree within 1 %. Possibly a revision "
                                             "or the same driver: not filled, not imported."})
                    summary.prefix_review += 1
                    summary.skipped.append((ref.hint("canonical_id"),
                                            f"prefix match for review: {old.get('driver_id')} (suffix {suffix!r})"))
                    continue
            if hit is not None:
                path, old = hit
                summary.matched += 1
                rows = compare_records(old, record)
                for row in rows:
                    conflicts.append({"driver_id": old.get("driver_id"), "path": str(path),
                                      "driver_base_url": record["field_provenance"]["fs_hz"]["url"]
                                      if "fs_hz" in record["field_provenance"] else ref.url,
                                      "note": _conflict_note(old, row["field"]), **row})
                summary.conflicts += len(rows)
                filled: List[str] = []
                if (fill and not result.incomplete and old.get("record_schema_version", 1) == 1
                        and path.resolve().is_relative_to(db_dir.resolve())):
                    updated, filled = fill_missing(old, record)
                    if updated is not None:
                        try:
                            _driver_from_dict(updated)
                        except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
                            filled = []
                            summary.skipped.append((ref.hint("canonical_id"), f"fill rejected by loader: {error}"))
                        else:
                            path.write_text(json.dumps(updated, indent=4) + "\n")
                            existing[key] = (path, updated)
                            summary.filled += 1
                            summary.outputs.append(str(path))
                matches.append({"canonical_id": ref.hint("canonical_id"), "driver_id": old.get("driver_id"),
                                "path": str(path), "matched_by": how,
                                "conflicting_fields": [r["field"] for r in rows], "filled": filled})
                continue
            if key in imported:
                summary.skipped.append((ref.hint("canonical_id"), f"duplicate of {imported[key]} in the dataset"))
                continue
            target = db_dir / record["manufacturer"] / f"{record['driver_id']}.json"
            if target.resolve() in taken_paths:
                impedance = (record["identity"]["variant"] or {}).get("impedance_ohm")
                record["driver_id"] = driver_id_for(record["manufacturer"], record["model_name"], impedance)
                target = db_dir / record["manufacturer"] / f"{record['driver_id']}.json"
                if target.resolve() in taken_paths:
                    summary.skipped.append((ref.hint("canonical_id"), f"driver_id {record['driver_id']} already used"))
                    continue
            imported[key] = ref.hint("canonical_id")
            reason = result.reason if result.incomplete else None
            if reason is None:
                errors, _ = check_structure(record)
                if errors:
                    reason = "schema: " + "; ".join(errors)
            if reason is None:
                try:
                    _driver_from_dict(record)
                except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError) as error:
                    reason = f"loader rejects record: {error}"
            if reason is None:
                path = write_record(db_dir, record)
                taken_paths.add(path.resolve())
                summary.outputs.append(str(path))
                summary.written += 1
            elif result.incomplete or record.get("driver_type") == "compression":
                summary.outputs.append(str(write_incomplete(incomplete_dir, record, reason)))
                summary.incomplete += 1
            else:
                summary.skipped.append((ref.hint("canonical_id"), reason))
                cls = reason.split(":")[0]
                summary.skip_reasons[cls] = summary.skip_reasons.get(cls, 0) + 1
        if summary.pages:
            write_reports(report_dir, summary, conflicts, matches, review)
        return summary


def _conflict_note(existing: dict, field: str) -> str:
    if field == "mms_kg" and any("Mmd" in str(n) for n in (existing.get("provenance") or {}).get("known_limitations", [])):
        return "existing mms_kg may be the dry mass Mmd"
    return ""


def write_reports(report_dir: Path, summary: ImportSummary, conflicts: List[dict], matches: List[dict],
                  review: Sequence[dict] = ()) -> None:
    """conflicts.csv, prefix_review.csv, matches.json and summary.json under *report_dir*."""
    report_dir.mkdir(parents=True, exist_ok=True)
    with (report_dir / "prefix_review.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREFIX_REVIEW_COLUMNS)
        writer.writeheader()
        for row in sorted(review, key=lambda r: (r["path"], r["canonical_id"])):
            writer.writerow({k: row.get(k) for k in PREFIX_REVIEW_COLUMNS})
    with (report_dir / "conflicts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONFLICT_COLUMNS)
        writer.writeheader()
        for row in sorted(conflicts, key=lambda r: (r["path"], r["field"])):
            writer.writerow({k: row.get(k) for k in CONFLICT_COLUMNS})
    (report_dir / "matches.json").write_text(json.dumps(sorted(matches, key=lambda m: m["path"]), indent=1) + "\n")
    reasons: Dict[str, int] = {}
    for _, reason in summary.skipped:
        cls = reason.split(":")[0]
        reasons[cls] = reasons.get(cls, 0) + 1
    (report_dir / "summary.json").write_text(json.dumps({
        "dataset": DATASET_URL, "license": LICENSE, "attribution": ATTRIBUTION,
        **{k: v for k, v in summary.to_dict().items() if k not in ("outputs", "skip_reasons")},
        "skipped_by_reason": dict(sorted(reasons.items(), key=lambda item: -item[1])),
    }, indent=1) + "\n")


__all__ = ["DriverBaseSource", "convert_entry", "compare_records", "consistency_problems", "fill_missing",
           "same_driver_by_values", "prefix_candidate", "DATASET_URL", "ATTRIBUTION"]
