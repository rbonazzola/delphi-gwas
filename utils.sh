#!/bin/bash

inverse_map() {
  local global=$1

  # bloque = parte entera de (global / 120)
  local block=$(( global / 120 ))
  local embedding=$(( global % 120 ))

  case $block in
    0) age=20 ;;
    1) age=30 ;;
    2) age=40 ;;
    3) age=50 ;;
    4) age=60 ;;
    *) echo "Índice inválido" >&2; return 1 ;;
  esac

  echo "$embedding $age"
}

