#!/bin/bash

# Usage:
#   ./Mdelta.bash "[1,2,3]" "[1,0,0]" mat01.txt

# Change this to the directory the Octate scripts live in
OCTAVE_SCRIPT_DIR="/home/paul/lcn/20251008-autoregoster-generate-offset-mats"

OCTAVE_FN="compute_MdeltaPrescription_fromSS_v01_RAS"

POS=$1
SS_dir=$2
OUTFILE=$3

if [ -z "$POS" ] || [ -z "$SS_dir" ]; then
  echo "Usage: $0 POS SS_dir [OUTFILE]"
  echo ""
  echo "Arguments:"
  echo "  POS     - Position vector in RAS coordinates, e.g., \"[1,2,3]\""
  echo "  SS_dir  - Slice-select direction in RAS coordinates, e.g., \"[1,0,0]\""
  echo "  OUTFILE - (optional) Output filename, e.g., \"mat01.txt\""
  echo ""
  echo "Example:"
  echo "  $0 \"[1,2,3]\" \"[1,0,0]\" mat01.txt"
  echo "  $0 \"[1,2,3]\" \"[1,0,0]\"  # outputs to stdout"
  exit 1
fi

if [ -z "$OUTFILE" ]; then
  EVAL_STR="compute_MdeltaPrescription_fromSS_v01_RAS(${POS}, ${SS_dir}); exit"
else
  EVAL_STR="compute_MdeltaPrescription_fromSS_v01_RAS(${POS}, ${SS_dir}, '$OUTFILE'); exit"
fi

echo "   POS is (RAS): $POS"
echo "SS_dir is (RAS): $SS_dir"
echo "     Outfile is: $OUTFILE"
#echo "    Eval str is: $EVAL_STR"

octave -p $OCTAVE_SCRIPT_DIR --eval "$EVAL_STR"
