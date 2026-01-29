#!/bin/bash
# gdb対話型デバッグスクリプト

cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf

# gdbで対話的に実行
gdb ./spa_mpi << 'GDB_EOF'
set args < in.surf.rotate
break Particle::sort
run
# クラッシュ時に自動的に停止
# 変数の値を確認
print grid->maxlocal
print grid->nlocal
print grid->nghost
print i
print particles[i].icell
print fixed_maxlocal
print fixed_nlocal
print fixed_nglocal
bt
info registers
continue
GDB_EOF
