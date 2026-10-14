"""Build the arXiv upload bundle for the paper (preprint.tex).

    cd paper && python make_arxiv.py

Writes paper/arxiv/ (a flat directory) and paper/arxiv.tar.gz, which is what gets
uploaded. arXiv does not run BibTeX, so the .bbl is generated here and shipped; it
must carry the same basename as the .tex. arXiv also publishes the source files,
so full-line and trailing comments are stripped from the .tex first: the internal
notes about venues and the ledger are for this repo, not for readers.

Re-run after any edit to preprint.tex, refs.bib or figures/.
"""
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "arxiv"
MAIN = "preprint"
SUPPORT = ["SciPost.cls", "SciPost_bibstyle.bst"]


def strip_comments(tex: str) -> str:
    out = []
    for line in tex.splitlines():
        if line.lstrip().startswith("%"):
            continue
        # Trailing comment preceded by whitespace and not escaped: "foo   % note" -> "foo".
        # A bare trailing "%" (e.g. "\textbf{%") suppresses a space and is kept.
        line = re.sub(r"(?<!\\)\s+%\s+\S.*$", "", line)
        out.append(line)
    text = "\n".join(out) + "\n"
    return re.sub(r"\n{3,}", "\n\n", text)


def main():
    # Fresh build so the .bbl matches the current refs.bib.
    subprocess.run(["tectonic", "-X", "compile", "--keep-intermediates", f"{MAIN}.tex"],
                   cwd=HERE, check=True)

    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "figures").mkdir(parents=True)

    tex = (HERE / f"{MAIN}.tex").read_text()
    (OUT / f"{MAIN}.tex").write_text(strip_comments(tex))
    shutil.copy(HERE / f"{MAIN}.bbl", OUT)
    for f in SUPPORT:
        shutil.copy(HERE / f, OUT)

    figs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
    for fig in figs:
        shutil.copy(HERE / fig, OUT / fig)

    tar = HERE / "arxiv.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        for p in sorted(OUT.rglob("*")):
            if p.is_file():
                t.add(p, arcname=str(p.relative_to(OUT)))

    for ext in ("aux", "log", "blg", "out", "toc", "bbl"):
        (HERE / f"{MAIN}.{ext}").unlink(missing_ok=True)
    print(f"{len(figs)} figures; wrote {OUT.relative_to(HERE)}/ and {tar.name}")


if __name__ == "__main__":
    main()
