#!/bin/bash

# script to regenerate the python documentation for all the source code in ../src.
# please run from the 'doc' directory

mkdir -p pydoc
rm pydoc/*.py.txt

cd ../src
for py in *.py; do
    pydoc ./$py > ../doc/pydoc/$py.txt
done
