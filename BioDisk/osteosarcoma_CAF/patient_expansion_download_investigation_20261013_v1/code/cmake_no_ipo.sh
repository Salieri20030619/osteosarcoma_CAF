#!/bin/sh
case "$1" in --build|--version|--install|-E) exec /usr/bin/cmake "$@" ;; *) exec /usr/bin/cmake "$@" -DTBB_ENABLE_IPO=OFF ;; esac
