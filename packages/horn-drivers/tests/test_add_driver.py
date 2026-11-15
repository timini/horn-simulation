"""horn-add-driver: datasheet values with explicit units in, validated SI record out."""
import json

import pytest

from horn_drivers.add_driver import main
from horn_drivers.loader import load_driver, load_drivers

SOURCE = "https://example.org/acme-ac10.pdf"
BASE = ["--manufacturer", "Acme", "--model", "AC10 8ohm", "--type", "cone",
        "--fs-hz", "45", "--re-ohm", "5.6", "--bl-tm", "15.2", "--sd-cm2", "346",
        "--le-mh", "0.9", "--qms", "5.1", "--qes", "0.38", "--qts", "0.354",
        "--xmax-mm", "4", "--power-w", "300", "--f-low-hz", "50", "--f-high-hz", "3000",
        "--source-url", SOURCE, "--retrieved-on", "2026-11-15"]


def test_datasheet_values_are_converted_and_written_with_provenance(tmp_path):
    assert main(["--db", str(tmp_path), *BASE, "--vas-l", "52"]) == 0
    path = tmp_path / "Acme" / "acme-ac10-8ohm.json"
    record = json.loads(path.read_text())
    assert record["parameter_source"] == SOURCE
    assert record["provenance"]["retrieved_on"] == "2026-11-15"
    assert "vas_l" in record["provenance"]["supplied_fields"]
    assert record["provenance"]["derived_fields"] == ["mms_kg"]
    assert record["record_schema_version"] == 1
    driver = load_driver(str(tmp_path), "acme-ac10-8ohm")
    assert driver.sd_m2 == pytest.approx(0.0346)
    assert driver.le_h == pytest.approx(0.0009)
    assert driver.xmax_m == pytest.approx(0.004)
    assert driver.mms_kg == pytest.approx(0.04079, rel=2e-3)
    assert (driver.usable_f_low_hz, driver.usable_f_high_hz) == (50, 3000)
    assert driver.power_w == 300


def test_compression_driver_exit_diameter_becomes_throat_area(tmp_path):
    args = ["--db", str(tmp_path), "--manufacturer", "Acme", "--model", "CD1", "--type", "compression",
            "--fs-hz", "784", "--re-ohm", "6.09", "--bl-tm", "2.06", "--sd-cm2", "16.2", "--mms-g", "0.3",
            "--le-mh", "0.04", "--qms", "2.02", "--qes", "2.03", "--exit-diameter-mm", "25.4",
            "--source-url", SOURCE]
    assert main(args) == 0
    (driver,) = load_drivers(str(tmp_path))
    assert driver.effective_throat_area == pytest.approx(3.14159265 * 0.0127 ** 2)


@pytest.mark.parametrize("extra,message", [
    (["--qts", "0.6"], "Q factor identity violated"),       # contradicts Qms and Qes
    (["--mms-g", "40000"], "mms_kg"),              # kilograms typed as grams
])
def test_inconsistent_or_mis_scaled_values_write_nothing(tmp_path, capsys, extra, message):
    args = ["--db", str(tmp_path), *BASE, *extra]
    if "--mms-g" not in extra:
        args += ["--vas-l", "52"]
    assert main(args) == 1
    assert message in capsys.readouterr().err
    assert not any(tmp_path.rglob("*.json"))


def test_source_url_is_required():
    with pytest.raises(SystemExit):
        main([a for a in BASE if a not in ("--source-url", SOURCE)])
