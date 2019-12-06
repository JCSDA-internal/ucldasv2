#!/bin/bash
set -e

if [[ "$#" -ne 1 ]]; then
  echo "Usage: ./rename_project.sh <new_proj_name>"
  exit 1
fi
name_new=$1
name_orig=mymodel

# make sure in lower case
name_new=${name_new,,}

# rename directories/files
mv src/$name_orig src/$name_new
mv cmake/mymodel_compiler_flags.cmake cmake/${name_new}_compiler_flags.cmake

# replace "mymodel" with correct name in files, respecting the case
name_orig_uc=${name_orig^^}
name_new_uc=${name_new^^}
files=$(grep $name_orig * -Rl --ignore-case)
for f in $files; do
  [[ "$f" == bundle/*/* ]] && continue
  [[ "$f" == "rename_project.sh" ]] && continue
  echo "replacing \"$name_orig\" in $f"  
  sed "s/$name_orig/$name_new/g; s/$name_orig_uc/$name_new_uc/g" $f > $f.tmp && mv $f.tmp $f
done