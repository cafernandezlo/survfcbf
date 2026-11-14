#!/bin/bash
# survfcbf Application Note: compile the three PDFs in paper/appnote/
#   appnote.pdf        the note in the journal's format (OUP authoring template)
#   supplementary.pdf  the supplementary material, one file for the journal and bioRxiv
#   preprint.pdf       the preprint: same text in a neutral format with line numbers
#
# Usage:
#   cd paper/appnote && bash build_appnote.sh            # all three
#   cd paper/appnote && bash build_appnote.sh appnote    # just one
#
# The tables of the supplementary material come from generated/, written by
#   Rscript --no-init-file paper/make_supplementary_tables.R
#
# Official OUP template (oup-authoring-template, CTAN). On a minimal TinyTeX the
# dependencies install once with:
#   tlmgr install oup-authoring-template crop sttools totcount arydshln subfloat \
#         jknapltx rsfs algorithms algorithmicx silence anyfontsize xurl \
#         newtx orcidlink lineno setspace
set -euo pipefail

TEXBIN="${TEXBIN:-$HOME/AppData/Roaming/TinyTeX/bin/windows}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

build() {
  local doc="$1"
  "$TEXBIN/pdflatex" -interaction=nonstopmode "$doc.tex" >/dev/null || true
  "$TEXBIN/bibtex" "$doc" >/dev/null
  "$TEXBIN/pdflatex" -interaction=nonstopmode "$doc.tex" >/dev/null || true
  # repeat until the contents and cross-references settle (a stale .toc needs more passes)
  local pass=0
  while [ $pass -lt 4 ]; do
    "$TEXBIN/pdflatex" -interaction=nonstopmode "$doc.tex" >/dev/null || true
    pass=$((pass + 1))
    grep -q "Rerun to get\|Label(s) may have changed" "$doc.log" || break
  done
  grep -h "^Output written" "$doc.log" || true

  # fail on a LaTeX error, or on unresolved references or citations
  if grep -q "^! " "$doc.log"; then
    grep -A3 "^! " "$doc.log" >&2
    echo "ERROR: $doc.tex does not compile" >&2
    exit 1
  fi
  if grep -E "Warning: (Citation|Reference) .* undefined" "$doc.log"; then
    echo "ERROR: unresolved references in $doc" >&2
    exit 1
  fi
  echo "OK: $doc.pdf"
}

for doc in "${@:-appnote supplementary preprint}"; do
  for d in $doc; do build "$d"; done
done

# the PDFs are fine for review with placeholder DOIs, never for submission or bioRxiv
if grep -qE "zenodo\.(XXXXXXX|YYYYYYY)" macros.tex; then
  echo "WARNING: macros.tex still holds placeholder Zenodo DOIs; do not submit or post these PDFs" >&2
fi
