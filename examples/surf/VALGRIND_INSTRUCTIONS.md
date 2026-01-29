# Valgrind調査の手順

## 1. Valgrindのインストール

Valgrindがインストールされていない場合は、以下のコマンドでインストールしてください:

```bash
sudo apt install valgrind
```

## 2. Valgrind調査の実行

インストール後、以下のいずれかの方法でValgrind調査を実行できます:

### 方法1: スクリプトを使用（推奨）

```bash
cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf
./run_valgrind.sh
```

### 方法2: 直接コマンドを実行

```bash
cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf
valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose --error-limit=no --log-file=valgrind.log ./spa_mpi < in.surf.rotate
```

## 3. 結果の確認

Valgrind調査が完了すると、`valgrind.log`ファイルに結果が保存されます。

主要なエラーを確認するには:

```bash
grep -A 5 'Invalid read\|Invalid write\|Use of uninitialised' valgrind.log
```

セグメンテーションフォルトの原因を特定するには:

```bash
grep -A 10 'SIGSEGV' valgrind.log
```

## 4. 現在の問題

現在、`Particle::sort()`内でセグメンテーションフォルトが発生しています。
- アドレス: +0x1404a6
- 問題: `cinfo[icell].first`へのアクセス時に無効なアドレスにアクセス
- `icell = 380`が`nglocal`を超えている可能性

Valgrind調査により、メモリ破壊や範囲外アクセスの詳細な情報が得られるはずです。
