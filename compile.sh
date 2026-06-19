#!/bin/bash
set -e
TEX_FILE="tex/ProjectedYukawaOperator.tex"; OUTPUT_DIR="out"; MAIN_NAME="ProjectedYukawaOperator"
mkdir -p "$OUTPUT_DIR"
export TEXINPUTS=".:./tex:${TEXINPUTS}"
pdflatex -file-line-error -interaction=nonstopmode -output-directory="$OUTPUT_DIR" "$TEX_FILE"
( cd "$OUTPUT_DIR" && BSTINPUTS="../tex:${BSTINPUTS}" BIBINPUTS="../tex:${BIBINPUTS}" bibtex "$MAIN_NAME" ) || true
pdflatex -file-line-error -interaction=nonstopmode -output-directory="$OUTPUT_DIR" "$TEX_FILE"
pdflatex -file-line-error -interaction=nonstopmode -output-directory="$OUTPUT_DIR" "$TEX_FILE"
echo "Done: $OUTPUT_DIR/$MAIN_NAME.pdf"
