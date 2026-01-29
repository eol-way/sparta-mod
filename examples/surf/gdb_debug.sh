#!/bin/bash
# gdbデバッグスクリプト

cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf

# gdbで実行
gdb -batch -ex "set args < in.surf.rotate" \
    -ex "run" \
    -ex "bt" \
    -ex "info registers" \
    -ex "print grid->maxlocal" \
    -ex "print grid->nlocal" \
    -ex "print grid->nghost" \
    -ex "print particles[i].icell" \
    -ex "print i" \
    -ex "print fixed_maxlocal" \
    -ex "print fixed_nlocal" \
    ./spa_mpi 2>&1 | tee gdb_output.log
