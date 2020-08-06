#!/bin/bash
CDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $CDIR/..

if [[ "$#" -ne 1 ]]; then
    echo "ERROR in usage: "
    echo " enable_all.sh < -1 | 0 | 1 >"
    echo "  -1 = undo the uncommenting"
    echo "   0 = normal usage, remove all \"TODO(template_impl)\" comments"
    echo "   1 = similar to 1, but comment out the comments, making it reversible"
    exit 1
fi

token="TODO(template_impl)"
files=$(grep "$token" {src,test}/ -Rl )

for f in $files; do
    echo $f
    if [[ "$1" == 1 ]]; then
        sed "s/\(#\\|\\/\\/\) \($token\) \(.*\)/\3 \1 \2/g" $f -i
        sed "s/\(\\/\\* $token\)/\\/\\/ \1/g" $f -i
        sed "s/\($token \\*\\/\)/\\/\\/ \1/g" $f -i

    elif [[ "$1" == 0 ]]; then
        sed "s/\(#\\|\\/\\/\) \($token\) \(.*\)/\3/g" $f -i
        sed "/\(\\/\\* $token\)/d" $f -i
        sed "/\($token \\*\\/\)/d" $f -i

    elif [[ "$1" == -1 ]]; then
        sed "s/\\/\\/ \(\\/\\* $token\)/\1/g" $f -i
        sed "s/\\/\\/ \($token \\*\\/\)/\1/g" $f -i
        sed "s/\([[:space:]]*\)\(.*\)\(#\\|\\/\\/\) \($token\)/\1\3 \4 \2/g" $f -i

    fi
done
