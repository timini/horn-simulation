"""Documentation checks: relative links resolve, and USAGE.md covers main.nf and the justfile."""
import re
import subprocess
from pathlib import Path
from urllib.parse import unquote

import pytest

ROOT = Path(__file__).resolve().parents[2]
USAGE = ROOT / "docs" / "USAGE.md"
SKIP_DIRS = {".git", ".venv", ".venv-cross-software", "research", "results", "work", ".nextflow", ".nf-test",
             ".tools", "node_modules"}
LINK = re.compile(r"!?\[[^\]]*\]\((<[^>]+>|[^)\s]+)(?:\s+\"[^\"]*\")?\)")


def markdown_files():
    """Tracked Markdown files, so local untracked notes do not fail the checks.

    Outside a git checkout, every Markdown file not under SKIP_DIRS.
    """
    try:
        listed = subprocess.run(["git", "ls-files", "-z", "*.md"], cwd=ROOT, capture_output=True,
                                text=True, check=True).stdout.split("\0")
        return sorted(ROOT / name for name in listed if name and (ROOT / name).is_file())
    except (OSError, subprocess.CalledProcessError):
        return [path for path in sorted(ROOT.rglob("*.md"))
                if not SKIP_DIRS.intersection(path.relative_to(ROOT).parts)]


def relative_links(text, with_anchor=False):
    """Relative link targets in Markdown text, outside fenced code blocks.

    Yields the path, or (path, anchor) with ``with_anchor``; a same-page link
    has an empty path.
    """
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    for match in LINK.finditer(text):
        target = match.group(1)
        if target.startswith("<"):
            target = target[1:-1]
        if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*:", target):
            continue
        path, _, anchor = target.partition("#")
        if with_anchor:
            yield path, anchor
        elif path:
            yield path


def heading_anchors(text):
    """GitHub-style anchors for the headings of a Markdown file."""
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    anchors, seen = set(), {}
    for heading in re.findall(r"^#{1,6}\s+(.+?)\s*#*\s*$", text, flags=re.M):
        slug = re.sub(r"[^\w\- ]", "", heading.replace("`", "").lower()).replace(" ", "-")
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        anchors.add(slug if count == 0 else f"{slug}-{count}")
    return anchors


def test_markdown_relative_links_resolve():
    broken = []
    files = markdown_files()
    assert any(p.name == "README.md" and p.parent == ROOT for p in files)
    for md in files:
        for target in relative_links(md.read_text(encoding="utf-8")):
            if not (md.parent / unquote(target)).exists():
                broken.append(f"{md.relative_to(ROOT)}: {target}")
    assert not broken, "Broken relative links:\n" + "\n".join(broken)


def test_markdown_anchors_resolve():
    broken = []
    for md in markdown_files():
        for path, anchor in relative_links(md.read_text(encoding="utf-8"), with_anchor=True):
            target = md if not path else md.parent / unquote(path)
            if not anchor or target.suffix != ".md" or not target.is_file():
                continue
            if anchor not in heading_anchors(target.read_text(encoding="utf-8")):
                broken.append(f"{md.relative_to(ROOT)}: {path}#{anchor}")
    assert not broken, "Links to missing headings:\n" + "\n".join(broken)


def test_heading_anchors_follow_github_rules():
    text = "# Run a study\n## Executed 800 to 1600 Hz example\n## `just run` (single)\n## Run a study"
    assert heading_anchors(text) == {"run-a-study", "executed-800-to-1600-hz-example", "just-run-single", "run-a-study-1"}


def test_link_parser_handles_angle_brackets_and_code():
    text = "[a](<docs/a b.md>) [b](c.md#x) [c](https://x.org) [d](#top)\n```\n[e](missing.md)\n```"
    assert list(relative_links(text)) == ["docs/a b.md", "c.md"]


def main_nf_params():
    """Top-level params.name = value assignments in main.nf, with the literal default."""
    params = {}
    for line in (ROOT / "main.nf").read_text().splitlines():
        match = re.match(r"^params\.(\w+)\s*=\s*(.+?)\s*(//.*)?$", line)
        if match:
            params[match.group(1)] = match.group(2).strip().strip('"').strip("'")
    return params


def usage_rows():
    rows = {}
    for line in USAGE.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\|\s*`(\w+)`\s*\|([^|]*)\|", line)
        if match:
            rows.setdefault(match.group(1), match.group(2).strip())
    return rows


def test_main_nf_params_found():
    params = main_nf_params()
    assert {"mode", "target_f_low", "outdir", "max_mouth_radius"} <= set(params)
    assert params["mode"] == "single"


@pytest.mark.parametrize("name,default", sorted(main_nf_params().items()))
def test_usage_documents_every_main_nf_param(name, default):
    rows = usage_rows()
    assert name in rows, f"docs/USAGE.md has no parameter row for `{name}`"
    cell = rows[name]
    if default == "null":
        assert cell.startswith("unset"), f"`{name}` defaults to null in main.nf; USAGE.md says {cell!r}"
    else:
        assert f"`{default}`" in cell, f"`{name}` defaults to {default} in main.nf; USAGE.md says {cell!r}"


def test_usage_lists_every_just_recipe():
    recipes = re.findall(r"^([a-z][\w-]*)(?:\s[^:]*)?:(?!=)", (ROOT / "justfile").read_text(), flags=re.M)
    assert "run-auto" in recipes and "test-local" in recipes
    text = USAGE.read_text(encoding="utf-8")
    missing = [r for r in recipes if r != "default" and f"`just {r}" not in text]
    assert not missing, f"docs/USAGE.md does not list: {missing}"
