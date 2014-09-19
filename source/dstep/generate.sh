#!/usr/bin/env bash

# requires installed dstep executable

STDARG_H_DIR=/usr/lib/clang/3.4.2/include # set this path to point to stdarg.h

clang -Xclang -E -fsyntax-only -DSAMTOOLS -I../c/htslib -I${STDARG_H_DIR} ../c/htslib/cram/cram.h > cram_all.h
rdmd remove_sysincludes.d cram_all.h
cat includes.h cram_all.h.filtered > cram_all.h
dstep -DSAMTOOLS -I../c/htslib -I${STDARG_H_DIR} cram_all.h  -o cram.d

