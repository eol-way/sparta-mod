# SPARTA ビルド環境 - Cursor引き継ぎ資料

## 【プロジェクト概要】

- **対象ソフト:** SPARTA (DSMCシミュレーションソフト)
- **環境:** Windows 11 + WSL2 (Ubuntu, Linux 6.6.87.2)
- **作業ディレクトリ:** `/home/yu0126016/SPARTA/GIT/sparta/`
- **目的:** 
  1. ソースコードの改造（物理モデルの追加等）
  2. ビルド・実行の自動化
  3. 結果の検証までをCursorで完結させる「バイブコーディング」の実践

## 【現在のステータス】✅

### 解決済みの問題
- ✅ **ビルドエラー解決:** `make linux`でエラーが出ていたが、WSL環境では`mpi`ターゲットを使用することで解決
- ✅ **compile_commands.json生成:** `bear -- make mpi`で正常に生成済み（`src/compile_commands.json`）
- ✅ **実行ファイル生成:** `spa_mpi`が正常にビルド済み（`src/spa_mpi`）

### 重要な発見
- **`linux`ターゲットは存在しない:** `MAKE/Makefile.linux`が存在しないため、`make linux`は失敗する
- **WSL環境での正しいターゲット:** `mpi`（MPI対応）または`serial`（シリアルビルド）
- **MPI環境:** `/usr/bin/mpicc`と`/usr/bin/mpicxx`が利用可能

## 【ビルドコマンド】

### 基本ビルド（compile_commands.json生成付き）
```bash
cd /home/yu0126016/SPARTA/GIT/sparta/src
bear -- make mpi
```

### クリーンビルド
```bash
cd /home/yu0126016/SPARTA/GIT/sparta/src
make clean-all
bear -- make mpi
```

### 利用可能なターゲット
- `mpi` - MPI対応ビルド（推奨、WSL環境で使用中）
- `serial` - シリアルビルド（MPIなし）
- `openmpi` - OpenMPI使用
- `mpi_big` - 大きな問題用のMPIビルド
- その他多数（`ls src/MAKE/Makefile.*`で確認可能）

## 【ディレクトリ構造】

```
/home/yu0126016/SPARTA/GIT/sparta/
├── src/                    # ソースコード
│   ├── Makefile           # メインMakefile
│   ├── compile_commands.json  # clangd用（bearで生成）
│   ├── spa_mpi            # ビルド済み実行ファイル
│   ├── Obj_mpi/           # ビルドオブジェクトファイル
│   ├── MAKE/              # マシン別Makefile
│   │   ├── Makefile.mpi   # MPI用設定（使用中）
│   │   ├── Makefile.serial
│   │   └── ...
│   └── *.cpp, *.h         # ソースファイル
└── build/                 # CMakeビルドディレクトリ（別用途）
```

## 【Cursorでの作業フロー】

### 1. コード解析環境の確認
```bash
# compile_commands.jsonが存在するか確認
ls -la /home/yu0126016/SPARTA/GIT/sparta/src/compile_commands.json
```

### 2. コード改造の実施
- CursorのComposer (`Ctrl + I`) を使用
- 例: 「`src/`以下のファイルを参考にして、新しい物理モデルを追加して」

### 3. ビルドとテスト
```bash
cd /home/yu0126016/SPARTA/GIT/sparta/src
bear -- make mpi
# 実行ファイルでテスト
./spa_mpi < input.in
```

### 4. エラー修正のサイクル
- ビルドエラー → Cursorで修正 → 再ビルド → テスト

## 【重要な注意点】

1. **ターゲット名:** WSL環境では`linux`ではなく`mpi`を使用
2. **compile_commands.json:** `bear`コマンドで生成（clangd用）
3. **ビルドディレクトリ:** `Obj_mpi/`にオブジェクトファイルが生成される
4. **実行ファイル:** `src/spa_mpi`が生成される

## 【トラブルシューティング】

### エラー: `make: *** [Makefile:96: linux] Error 1`
- **原因:** `MAKE/Makefile.linux`が存在しない
- **解決:** `make mpi`を使用

### compile_commands.jsonが生成されない
- `bear`がインストールされているか確認: `which bear`
- ビルドが成功しているか確認（エラーがあると生成されない場合がある）

### MPI関連エラー
- MPIがインストールされているか確認: `which mpicc mpicxx`
- インストール: `sudo apt-get install libopenmpi-dev openmpi-bin`

## 【次のステップ】

1. ✅ ビルド環境の構築（完了）
2. ✅ compile_commands.jsonの生成（完了）
3. ⏭️ コード改造の実施
4. ⏭️ テストの自動化
5. ⏭️ 結果の検証

---

**最終更新:** 2024年12月26日
**状態:** ビルド環境構築完了、開発準備完了




