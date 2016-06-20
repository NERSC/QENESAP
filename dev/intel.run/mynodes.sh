#!/bin/sh

NODES="localhost"
#NODES="$(bjobs | tail -n +2 | tr -s ' ' | cut -d' ' -f1-6 | rev | cut -d' ' -f1 | rev)"

INDEX=0
if [ "" != "$1" ]; then
  INDEX=$1
  shift
fi

COUNT=0
for IP in ${NODES}; do
  if [ "0" != "$((INDEX<=COUNT))" ]; then
    echo ${IP}
  fi
  COUNT=$((COUNT+1))
done

