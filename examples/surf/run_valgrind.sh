#!/bin/bash
# Valgrind調査用スクリプト
# 使用方法: ./run_valgrind.sh

cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf

# Valgrindがインストールされているか確認
if ! command -v valgrind &> /dev/null; then
    echo "エラー: Valgrindがインストールされていません"
    echo "以下のコマンドでインストールしてください:"
    echo "sudo apt install valgrind"
    exit 1
fi

# Valgrind調査を実行
echo "Valgrind調査を開始します..."
valgrind \
    --tool=memcheck \
    --leak-check=full \
    --show-leak-kinds=all \
    --track-origins=yes \
    --verbose \
    --error-limit=no \
    --log-file=valgrind.log \
    ./spa_mpi < in.surf.rotate

echo ""
echo "Valgrind調査が完了しました。結果は valgrind.log に保存されています。"
echo "主要なエラーを確認するには:"
echo "grep -A 5 'Invalid read\|Invalid write\|Use of uninitialised' valgrind.log"
