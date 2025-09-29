#! /bin/bash -x

cd `dirname $0`
script=1dmodel_analysis.py

### Visalization ###
echo "+visualize by python"

if type "python" > /dev/null 2>&1; then
    python ${script}
elif type "python3" > /dev/null 2>&1; then
    python3 ${script}
else
    echo "python does not found"
fi
