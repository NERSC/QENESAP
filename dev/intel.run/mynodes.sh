#!/bin/sh

#NODES="knl1 knl2 knl3 knl4 knl5 knl6 knl7 knl8"
NODES="localhost"

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

