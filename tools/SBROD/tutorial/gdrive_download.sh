#!/usr/bin/env bash

URL=$1
#ggID=$(echo $URL | sed 's/[^A-Za-z0-9_-]/\n/g' | sed -rn '/.{26}/p')
ggID=$(echo $URL | sed 's/=/ /g' | cut -d' ' -f2)
ggURL='https://drive.google.com/uc?export=download'
curl -sc /tmp/gcokie "${ggURL}&id=${ggID}" >/dev/null
getcode="$(awk '/_warning_/ {print $NF}' /tmp/gcokie)"
if [[ !  -z  $getcode  ]]; then code="&confirm="$getcode; else code=""; fi
curl -LOJb /tmp/gcokie "${ggURL}${code}&id=${ggID}"
