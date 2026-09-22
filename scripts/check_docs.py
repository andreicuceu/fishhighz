"""Check that the built API reference contains signatures and link targets."""

import argparse
import re
from html.parser import HTMLParser
from pathlib import Path


class APIPage(HTMLParser):
    """Collect object anchors, rendered signatures, and visible HTML text."""

    def __init__(self):
        super().__init__()
        self.anchors = set()
        self.signatures = set()
        self.text = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "id" in attributes:
            self.anchors.add(attributes["id"])
            if tag == "dt" and "sig" in attributes.get("class", "").split():
                self.signatures.add(attributes["id"])

    def handle_data(self, data):
        self.text.append(data)


def check_api(source, html):
    """Return errors for absent API targets or unparsed autodoc output."""
    errors = []
    total = 0
    for page in sorted(source.glob("*.md")):
        text = page.read_text(encoding="utf-8")
        # Recognize both spellings so this check also detects the original bug.
        targets = re.findall(
            r"(?:\.\. auto(?:class|function)::|```\{auto(?:class|function)\})"
            r"\s+(fishhighz\.[\w.]+)",
            text,
        )
        if not targets:
            continue
        path = html / "api" / page.with_suffix(".html").name
        if not path.is_file():
            errors.append(f"Missing API page: {path}")
            continue
        parsed = APIPage()
        parsed.feed(path.read_text(encoding="utf-8"))
        for target in targets:
            total += 1
            if target not in parsed.anchors or target not in parsed.signatures:
                errors.append(f"{path.name}: missing signature/anchor for {target}")
        for block in re.findall(r"```\{eval-rst\}\n(.*?)```", text, re.S):
            match = re.search(r"\.\. autoclass::\s+(\S+)", block)
            if match:
                for members in re.findall(r":(?:special-)?members:\s*([^\n]+)", block):
                    for member in members.split(","):
                        target = f"{match[1]}.{member.strip()}"
                        if target not in parsed.signatures:
                            errors.append(f"{path.name}: missing member {target}")
        if re.search(r"\.\.\s+py:\w+::", "".join(parsed.text)):
            errors.append(f"{path.name}: raw Python-domain directive in HTML")
    if not total:
        errors.append("No API targets checked")
    return total, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", type=Path, help="Sphinx HTML output directory")
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[1] / "docs" / "api"
    total, errors = check_api(source, args.html)
    if errors:
        print("\n".join(errors))
        raise SystemExit(1)
    print(f"Verified {total} API objects and their selected member signatures/anchors.")


if __name__ == "__main__":
    main()
