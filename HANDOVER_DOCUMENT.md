# SPARTA改造プロジェクト：引き継ぎ資料

**作成日**: 2024年12月26日  
**プロジェクト**: SPARTA移動壁・回転壁機能の実装  
**フェーズ**: Phase 1 月1完了（検知ロジック実装済み）

---

## 📋 目次

1. [プロジェクト概要](#プロジェクト概要)
2. [実施済み作業](#実施済み作業)
3. [実装内容の詳細](#実装内容の詳細)
4. [ファイル構成](#ファイル構成)
5. [ビルド環境](#ビルド環境)
6. [検証結果](#検証結果)
7. [次のステップ](#次のステップ)
8. [重要な注意事項](#重要な注意事項)
9. [参考資料](#参考資料)
10. [よくある質問（FAQ）](#よくある質問faq)

---

## プロジェクト概要

### 目的
DSMCコード「SPARTA」を改良し、回転式容積ポンプ解析に必要な「移動壁・回転壁」の物理モデルを確立する。

### 本年度のゴール
- 壁移動に伴う粒子損失（食われ）をゼロにする「粒子押し出し機能」の完全実装
- 回転壁における「法線方向速度成分」を考慮した反射モデルの実装

### 開発計画（3カ年）

#### 1年目（Phase 1）: 基盤技術（CPU）の確立
- **テーマ**: "Correctness"（正しさ）
- **目標**: 移動壁による質量損失をゼロにするアルゴリズムの実装と基礎検証
- **主要タスク**:
  - 月1: 検知ロジックの実装（✅ 完了）
  - 月2: 押し出しロジックの実装（⏭️ 未着手）
  - 月3: 並列化対応（MPI通信境界の処理）

#### 2年目（Phase 2）: 2D実機解析とGPU化・設計部門への移管
- **テーマ**: "Speed & Release"（速度とリリース）
- **目標**: 実機形状での解析検証とGPU（KOKKOS）への高速化実装
- **重要マイルストーン**: 年度末に「2D GPU版」を設計部門へ正式リリース

#### 3年目（Phase 3）: 3次元化と高度化
- **テーマ**: "Expansion"（3D拡張）
- **目標**: 3D解析および熱連成への拡張

### 現在のステータス
- ✅ **Phase 1 月1完了**: 粒子検知ロジックの実装と検証
- ⏭️ **Phase 1 月2**: 押し出しロジックの実装（未着手）

### 技術的採用手法（アルゴリズム方針）

#### 粒子押し出し処理（Phase 1: Wall Move & Sweep）
**手法**: 壁移動更新直後に、壁内部（Solid領域）に取り残された粒子を検知し、救済する。

**アルゴリズム**:
1. **Cell Selection**: 粒子が存在するセルの表面要素のみを対象とする
2. **Velocity Filter**: $V_n = \vec{v}_{wall} \cdot \vec{n} > 0$ （押し出し方向）の面のみ候補とする
3. **Projection Check**: 射影パラメータ $0 \le t \le 1$ により、線分範囲外の誤判定を防ぐ
4. **Re-positioning**: 最短距離 $|D|$ の面に、安全マージン $\epsilon$ を加えて再配置する
5. **Velocity Update**: 押し出し後にCLLモデル（または他の表面衝突モデル）を使用して速度を更新する
   - 移動壁を考慮した反射速度が計算される
   - 壁から離れる方向に速度が設定される
   - これにより、壁が粒子より速い場合でも再取り込まれない

#### 衝突判定アルゴリズムの選定根拠
- **離散衝突検出（DCD）**: 高速な粒子（能動的衝突）に使用。計算が極めて速い
- **連続衝突検出（CCD）**: 計算コストが高すぎるため採用しない
- **押し出し（Sweeping）**: DCDが見逃した「轢かれた粒子」を救済。低コストな距離計算で実現

**結論**: 「速い粒子はDCDで捌き、遅い粒子（壁に轢かれるやつ）は事後の押し出し処理で救う」という二段構えにより、CCD並みの質量保存精度を、DCD並みの軽さで実現する。

---

## 実施済み作業

### 1. 環境構築

#### ビルド環境の確立
- **WSL環境**: Linux 6.6.87.2-microsoft-standard-WSL2
- **ビルドコマンド**: `bear -- make mpi`（`src`ディレクトリで実行）
- **実行ファイル**: `src/spa_mpi`
- **compile_commands.json**: 生成済み（clangd用）

#### 重要な発見
- WSL環境では`make linux`ではなく`make mpi`を使用
- `MAKE/Makefile.linux`は存在しない
- MPI環境は利用可能（`/usr/bin/mpicc`, `/usr/bin/mpicxx`）

### 2. コード分析

#### SPARTAの構造理解
- **粒子の位置判定**: 2段階システム（セルレベル + 粒子レベル）
  - `Grid::set_inout()`: セルタイプ（OUTSIDE/INSIDE/OVERLAP）の分類
  - `Grid::outside_surfs()`: 個別粒子の壁内外判定
- **現在の粒子削除**: `move_surf::remove_particles()`で削除
- **表面衝突**: `surf_collide`クラスは既に移動壁速度計算機能を実装済み（反射速度計算用）

### 3. 計画の見直しと改善提案

#### SPARTAコード構造の分析結果
- **粒子の位置判定**: 2段階システム（セルレベル + 粒子レベル）
  - `Grid::set_inout()`: セルタイプ（OUTSIDE/INSIDE/OVERLAP/UNKNOWN）の分類
  - `Grid::outside_surfs()`: 個別粒子の壁内外判定
- **現在の粒子削除**: `move_surf::remove_particles()`でセル全体を削除対象とする（問題点あり）
- **表面衝突**: `surf_collide`クラスは既に移動壁速度計算機能を実装済み（反射速度計算用）

#### 改善提案（詳細は`PLAN_REVIEW.md`を参照）
1. **2段階検知システムの採用**: セルレベル + 粒子レベル（✅ 実装済み）
   - Stage 1: セルレベルで候補セルを絞り込み（高速化）
   - Stage 2: 粒子レベルで個別判定（精度確保）

2. **表面速度計算機能の追加**: 回転壁対応（⏭️ 月2で実装予定）
   - 回転壁: $\vec{v} = \vec{\omega} \times \vec{r} + \vec{v}_{trans}$
   - 各表面要素の速度を個別に計算する必要がある

3. **セル割り当て更新の実装**: 押し出し後の処理（⏭️ 月2で実装予定）
   - 押し出し後に`grid->assign_particles()`を呼ぶ
   - 必要に応じて`grid->acquire_ghosts()`も再実行

4. **動的安全マージンの使用**: グリッドサイズに応じた調整（⏭️ 月2で実装予定）
   - `grid->cell_epsilon`（最小セルサイズの半分）を使用
   - または`EPSSURF * min_edge_length`（表面要素の最小エッジ長）

#### 技術的リスクと対策
1. **凹部・エッジでの判定競合**: 速度フィルタ（$V_n > 0$）と射影パラメータ（$0 \le t \le 1$）の組み合わせで対処
2. **数値誤差による再削除**: 安全マージン $\epsilon$ を付与して再配置
3. **高速回転時のトンネリング**: $\Delta t$ の制約条件（ガイドライン）を明確化（本年度の範囲内）

### 4. 検知ロジックの実装

#### 実装した機能
- **関数**: `MoveSurf::detect_trapped_particles(int *&trapped_list, int &ntrapped)`
- **機能**: 2段階検知システム
  - Stage 1: セルレベルでINSIDEセルや移動した表面を含むセルを候補として取得
  - Stage 2: 各候補セル内の粒子を`grid->outside_surfs()`で個別判定
- **出力**: 検知された粒子数と最初の10個の粒子の詳細情報（ID、位置）

#### 統合箇所
- `fix_move_surf::end_of_step()`内で検知ロジックを呼び出し
- 削除処理の前に実行（現在は検知後に既存の削除処理が実行される）

### 5. 検証

#### テストケース
- **スクリプト**: `examples/surf/in.surf.slide`
- **内容**: 2つの円形表面が互いに接近するように移動
- **結果**: 正常動作確認（44〜68個の粒子が各ステップで検知される）

---

## 実装内容の詳細

### 修正したファイル

#### 1. `src/move_surf.h`
```cpp
// 追加した関数宣言
void detect_trapped_particles(int *&, int &);
```

#### 2. `src/move_surf.cpp`
```cpp
// 追加した関数実装
void MoveSurf::detect_trapped_particles(int *&trapped_list, int &ntrapped)
{
  // 2段階検知システム
  // 1. セルレベルで候補セルを取得
  // 2. 粒子レベルで個別判定
  // 3. ログ出力
}
```

**実装のポイント**:
- `grid->outside_surfs()`を使用して個別粒子を判定
- `grid->point_outside_surfs()`でセル内の参照点を取得
- MPI集約でグローバルな粒子数を計算
- メモリ動的確保（初期1000個、必要に応じて拡張）

#### 3. `src/fix_move_surf.cpp`
```cpp
// end_of_step()内に追加
int *trapped_list = NULL;
int ntrapped = 0;
if (particle->exist) {
  movesurf->detect_trapped_particles(trapped_list, ntrapped);
}
// ... 既存の削除処理 ...
if (trapped_list) memory->sfree(trapped_list);
```

**追加したinclude**:
- `#include "particle.h"`（`move_surf.cpp`に追加）

---

## ファイル構成

### プロジェクトルート
```
/home/yu0126016/SPARTA/GIT/sparta/
├── src/                          # ソースコード
│   ├── move_surf.h              # ✅ 修正済み
│   ├── move_surf.cpp            # ✅ 修正済み
│   ├── fix_move_surf.h          # 変更なし
│   ├── fix_move_surf.cpp        # ✅ 修正済み
│   ├── compile_commands.json    # bearで生成
│   ├── spa_mpi                  # ビルド済み実行ファイル
│   └── Obj_mpi/                 # ビルドオブジェクト
├── examples/
│   └── surf/
│       └── in.surf.slide        # 検証用スクリプト
└── [ドキュメントファイル]
    ├── PROJECT_PLAN.md          # 3カ年計画書
    ├── PLAN_REVIEW.md           # 計画レビューと改善提案
    ├── IMPLEMENTATION_PHASE1.md # Phase 1実装計画
    ├── DETECTION_IMPLEMENTATION.md  # 検知ロジック実装報告
    ├── DETECTION_VERIFICATION.md    # 検証結果
    ├── CURSOR_CONTEXT.md        # Cursor用コンテキスト（旧）
    └── HANDOVER_DOCUMENT.md     # このファイル
```

### 重要なディレクトリ
- **`src/`**: ソースコード（155個の.cpp、167個の.h）
- **`src/MAKE/`**: マシン別Makefile（`Makefile.mpi`を使用中）
- **`examples/`**: テストケース
- **`build/`**: CMakeビルドディレクトリ（別用途）

---

## ビルド環境

### ビルドコマンド
```bash
cd /home/yu0126016/SPARTA/GIT/sparta/src
bear -- make mpi
```

### 実行コマンド
```bash
cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf
/home/yu0126016/SPARTA/GIT/sparta/src/spa_mpi < in.surf.slide
```

### 環境情報
- **OS**: WSL2 (Ubuntu, Linux 6.6.87.2-microsoft-standard-WSL2)
- **コンパイラ**: mpic++ (MPI対応)
- **ビルドシステム**: Make（`src/Makefile`を使用）
- **ツール**: bear（compile_commands.json生成用）
- **MPI環境**: `/usr/bin/mpicc`, `/usr/bin/mpicxx`が利用可能

### 利用可能なビルドターゲット
- `mpi` - MPI対応ビルド（✅ 現在使用中、推奨）
- `serial` - シリアルビルド（MPIなし）
- `openmpi` - OpenMPI使用
- `mpi_big` - 大きな問題用のMPIビルド
- その他多数（`ls src/MAKE/Makefile.*`で確認可能）

### 注意事項
- **ターゲット名**: `mpi`を使用（`linux`は存在しない）
- **`MAKE/Makefile.linux`は存在しない**: WSL環境では`make linux`ではなく`make mpi`を使用
- **ワークスペースルート**: `/home/yu0126016/SPARTA/GIT/sparta/`に変更推奨
- **CMakeビルド**: 別途CMakeビルドシステムも利用可能（`BUILD_CMAKE.md`参照）が、現在はMakeを使用

### トラブルシューティング
- **エラー**: `make: *** [Makefile:96: linux] Error 1`
  - **原因**: `MAKE/Makefile.linux`が存在しない
  - **解決**: `make mpi`を使用
- **compile_commands.jsonが生成されない**
  - `bear`がインストールされているか確認: `which bear`
  - ビルドが成功しているか確認（エラーがあると生成されない場合がある）
- **MPI関連エラー**
  - MPIがインストールされているか確認: `which mpicc mpicxx`
  - インストール: `sudo apt-get install libopenmpi-dev openmpi-bin`

---

## 検証結果

### テストケース: `examples/surf/in.surf.slide`

#### 検証内容
- 2つの円形表面（group 1とgroup 2）が互いに接近するように移動
- `fix move/surf`コマンドで100ステップごとに表面を移動
- 各移動タイミングで検知ロジックが実行される

#### 結果
✅ **正常動作を確認**

**動作確認項目**:
- [x] 検知ロジックが実行される
- [x] 検知された粒子数が正しくカウントされる
- [x] 粒子の詳細情報（ID、位置）が正しく出力される
- [x] MPI集約が正常に動作する（1プロセスでも動作確認）
- [x] ログファイルとスクリーン出力の両方に出力される
- [x] 壁移動の各タイミングで検知が実行される

**検知結果の統計**:
- **検知メッセージ**: 各ステップで正常に出力
- **検知粒子数の範囲**: 37〜68個
  - 最小: 37個（step 1700）
  - 最大: 68個（step 1200）
  - 平均: 約52個（概算）
- **粒子詳細情報**: IDと位置座標が正しく出力
- **MPI集約**: 正常に動作

#### 出力例
```
Step 600:
Detected 56 trapped particles
  First 10 trapped particles (proc 0):
    Particle ID 1960394233: pos = (4.928255e+00, 5.981214e+00, 0.000000e+00)
    Particle ID 1295277939: pos = (4.941149e+00, 5.983245e+00, 0.000000e+00)
    ...

Step 700:
Detected 55 trapped particles
  First 10 trapped particles (proc 0):
    Particle ID 645513075: pos = (4.880441e+00, 5.927385e+00, 0.000000e+00)
    ...
```

#### 観察事項
1. **検知タイミング**: `fix move/surf`の`nevery=100`に合わせて、100ステップごとに検知が実行される
2. **粒子位置の妥当性**: 検知された粒子の位置は、移動している円形表面の近傍に集中している
3. **現在の動作**: 検知後、既存の`remove_particles()`が実行され、検知された粒子は削除される（Phase 1月1の目標達成）

詳細は`DETECTION_VERIFICATION.md`を参照。

---

## 次のステップ

### Phase 1 月2: 押し出しロジックの実装

#### タスク1: 表面速度計算機能の追加
**目標**: 回転・並進を考慮した表面速度計算機能を実装する

- **関数**: `MoveSurf::compute_surf_velocity(int isurf, double *vwall)`
- **機能**: 回転・並進を考慮した表面速度計算
  - 回転壁: $\vec{v} = \vec{\omega} \times \vec{r} + \vec{v}_{trans}$
  - 並進のみ: $\vec{v} = \vec{v}_{trans}$
- **実装箇所**: `src/move_surf.cpp`
- **参照**: `fix_move_surf`の引数から速度パラメータを取得

#### タスク2: 幾何計算関数の実装
**目標**: 点から表面要素への最短距離を計算する関数を実装する

- 点から線分への最短距離（2D）
- 点から三角形への最短距離（3D）
- 射影パラメータの計算（$0 \le t \le 1$）
- **注意**: 線分範囲外の誤判定を防ぐため、射影パラメータの範囲チェックが重要

#### タスク3: 押し出しアルゴリズムの実装
**目標**: 検知された粒子を壁外に押し出す処理を実装する

- **関数**: `MoveSurf::push_out_particles()`
- **手順**:
  1. **Cell Selection**: 検知された粒子が存在するセルの表面要素を取得
     - `cells[icell].csurfs`を使用
  2. **Velocity Filter**: 各表面要素について速度フィルタリング
     - `compute_surf_velocity()`で速度計算
     - 法線ベクトルは`surf->lines/tris[isurf].norm`を使用
     - $V_n = \vec{v}_{wall} \cdot \vec{n} > 0$ で判定（押し出し方向のみ）
  3. **Projection Check**: 射影パラメータの範囲チェック
     - 粒子位置から表面への最短距離を計算
     - 射影パラメータ $t$ が $[0,1]$ の範囲内か確認
     - 2D: 線分への射影、3D: 三角形への射影
  4. **Re-positioning**: 最短距離の面に安全マージンを加えて再配置
     - 最短距離の表面を選択
     - 安全マージン: `epsilon = max(EPSSURF * min_edge, grid->cell_epsilon)`
     - 粒子位置を更新: `x_new = x_surf + epsilon * norm`
  5. **Velocity Update**: CLLモデル（または他の表面衝突モデル）を使用して速度を更新
     - 表面要素の衝突モデルを取得: `surf->sc[tri->isc]` または `surf->sc[line->isc]`
     - 法線ベクトルを取得: `tri->norm` または `line->norm`
     - `surf->sc[isc]->collide(&particles[ip], dtremain, isurf, norm, -1, reaction)` を呼び出す
     - 移動壁を考慮した反射速度が自動的に計算される（CLLモデルの内部処理）
     - **重要な理由**: 壁が粒子より速い場合、速度を更新しないと再取り込まれる可能性がある

#### タスク4: セル割り当ての更新
**目標**: 押し出し後に粒子を新しいセルに割り当てる

- 押し出し後に`grid->assign_particles()`を呼ぶ
- 必要に応じて`grid->acquire_ghosts()`も再実行
- プロセス間移動がある場合は`comm->migrate_particles()`も考慮

#### タスク5: 削除処理の置き換え
**目標**: 削除処理を押し出し処理に置き換える

- `fix_move_surf::end_of_step()`内で`remove_particles()`の代わりに`push_out_particles()`を呼ぶ
- 検知ロジックとの統合（`detect_trapped_particles()`で検知した粒子を押し出す）

### Phase 1 月3: 並列化対応

#### タスク1: 粒子のプロセス間移動処理
- 押し出しにより粒子が別プロセスのセルに移動する可能性
- `comm->migrate_particles()`の呼び出しを検討

#### タスク2: Ghost cell更新の統合
- 押し出し後に`grid->acquire_ghosts()`を再実行
- 通信処理のデバッグ

#### タスク3: メモリリークチェック
- 動的メモリ確保・解放の確認
- Valgrind等を使用したメモリチェック

#### タスク4: ベンチマークテスト（BM-1）
- **テストケース**: 1次元ピストン（高速）
- **検証項目**: 粒子数が保存され、理論通りの衝撃波面が形成されるか確認

詳細は`IMPLEMENTATION_PHASE1.md`と`PLAN_REVIEW.md`を参照。

---

## 重要な注意事項

### コード構造

#### 現在の処理フロー（`fix_move_surf::end_of_step()`）
1. 表面移動 (`movesurf->move_lines/tris()`)
2. グリッド再構築 (`grid->surf2grid()`, `grid->set_inout()`)
3. **検知ロジック** (`movesurf->detect_trapped_particles()`) ← 追加済み
4. 粒子削除 (`movesurf->remove_particles()`) ← 月2で押し出しに置き換え予定
5. 通知 (`grid->notify_changed()`)

#### 重要な関数
- `Grid::outside_surfs(int icell, double *x, double *xcell)`: 粒子の壁内外判定
  - 戻り値: 1=壁外、0=壁内
- `Grid::point_outside_surfs(int icell, double *x)`: セル内の参照点取得
- `MoveSurf::remove_particles()`: 現在は削除処理（月2で押し出しに置き換え）

### 実装上の注意点

1. **メモリ管理**
   - `trapped_list`は動的確保（`memory->smalloc()`）
   - 使用後は`memory->sfree()`で解放

2. **MPI処理**
   - 検知ロジック内でMPI集約を使用
   - ログ出力は`me == 0`でのみ実行

3. **粒子のアクセス**
   - `particle->particles[i]`: 粒子データ
   - `cinfo[icell].first`: セル内の最初の粒子インデックス
   - `particle->next[i]`: 次の粒子へのリンク

4. **安全マージン**
   - 現在は未使用（検知のみ）
   - 月2で`EPSSURF`または`grid->cell_epsilon`を使用予定

### 既知の制約事項

1. **分散サーフェス非対応**
   - `surf->distributed`が1の場合はエラー
   - `fix_move_surf`のコンストラクタでチェック済み

2. **現在は削除処理が実行される**
   - 検知後、既存の`remove_particles()`が実行される
   - 月2で押し出し処理に置き換える予定

---

## 参考資料

### ドキュメントファイル
本プロジェクトで作成されたドキュメント一覧：

1. **`PROJECT_PLAN.md`**: 3カ年開発計画書（全体像、技術的採用手法、中長期ロードマップ）
2. **`PLAN_REVIEW.md`**: 計画レビューと改善提案（重要、技術的リスクの再評価）
3. **`IMPLEMENTATION_PHASE1.md`**: Phase 1の実装詳細計画（月1〜月3のタスク詳細）
4. **`DETECTION_IMPLEMENTATION.md`**: 検知ロジック実装報告（実装内容の詳細）
5. **`DETECTION_VERIFICATION.md`**: 検証結果の詳細（テストケースの詳細な結果）
6. **`BUILD_CMAKE.md`**: CMakeビルド情報（参考、現在はMakeを使用）
7. **`CURSOR_CONTEXT.md`**: Cursor用コンテキスト（旧版、このドキュメントに統合済み）
8. **`HANDOVER_DOCUMENT.md`**: このファイル（引き継ぎ資料）

### コード参照
重要なソースファイルと実装箇所：

- **`src/move_surf.h`**: `detect_trapped_particles()`関数の宣言
- **`src/move_surf.cpp`**: 検知ロジックの実装（約895行目以降、`detect_trapped_particles()`関数）
- **`src/fix_move_surf.cpp`**: Fixクラスのラッパー（134行目以降、`end_of_step()`関数）
- **`src/grid_surf.cpp`**: 粒子判定関数（`outside_surfs()`, `point_outside_surfs()`）

### テストケース
- **`examples/surf/in.surf.slide`**: 検証用スクリプト（2つの円形表面が接近するテスト）
- **`examples/surf/in.surf.rotate`**: 回転壁のテストケース（コメントアウト済み、将来使用予定）

### 3カ年計画の詳細
詳細な開発計画、技術的リスク、アルゴリズム選定根拠については`PROJECT_PLAN.md`を参照してください。
- Phase 1（1年目）: 基盤技術の確立
- Phase 2（2年目）: 2D GPU化と設計部門への移管
- Phase 3（3年目）: 3D拡張と高度化

---

## よくある質問（FAQ）

### Q1: ビルドエラーが出る
**A**: `make mpi`を使用してください。`make linux`は存在しません。WSL環境では`mpi`ターゲットを使用してください。

### Q2: 検知ロジックが実行されない
**A**: `fix move/surf`コマンドが正しく設定されているか確認してください。`nevery`パラメータで実行頻度を制御します。`end_of_step()`が呼ばれるタイミングで実行されます。

### Q3: 検知された粒子数が0になる
**A**: 壁移動が発生していない可能性があります。`fix move/surf`の設定を確認してください。また、壁が移動していない初期ステップでは検知されない可能性があります。

### Q4: ワークスペースでファイルが見つからない
**A**: ワークスペースのルートを`/home/yu0126016/SPARTA/GIT/sparta/`に設定してください。

### Q5: 検知ロジックの実装箇所はどこか
**A**: `src/move_surf.cpp`の約895行目以降に`detect_trapped_particles()`関数が実装されています。`src/move_surf.h`に関数宣言があります。

### Q6: 押し出しロジックの実装はいつ完了するか
**A**: Phase 1 月2で実装予定です。詳細は「次のステップ」セクションを参照してください。

### Q7: 現在の検知ロジックと削除処理の関係は
**A**: 現在は検知ロジックの後に既存の`remove_particles()`が実行され、検知された粒子は削除されます。Phase 1 月2で押し出し処理に置き換える予定です。

### Q8: 分散サーフェス（distributed surf）は対応しているか
**A**: 現在は非対応です。`surf->distributed`が1の場合はエラーになります。`fix_move_surf`のコンストラクタでチェック済みです。

---

## 連絡先・サポート

### プロジェクト情報
- **プロジェクト名**: SPARTA改造プロジェクト
- **目的**: 移動壁・回転壁機能の実装
- **現在のフェーズ**: Phase 1 月1完了

### 次の担当者へのメッセージ

1. **検知ロジックは正常に動作しています**
   - 実装済みで検証も完了
   - コードは`src/move_surf.cpp`に追加済み

2. **次のステップは押し出しロジックの実装**
   - `PLAN_REVIEW.md`の改善提案を参照
   - 表面速度計算が重要

3. **既存コードの理解**
   - `Grid::outside_surfs()`の動作を理解する
   - `fix_move_surf::end_of_step()`の処理フローを把握する

4. **質問があれば**
   - このドキュメントの関連セクションを参照
   - 実装済みコードを確認

---

**最終更新**: 2024年12月26日  
**ステータス**: Phase 1 月1完了、月2準備完了

