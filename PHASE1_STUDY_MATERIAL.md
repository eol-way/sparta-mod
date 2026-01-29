# Phase 1: 基礎理解 - 勉強資料

**作成日**: 2024年12月26日  
**対象**: Phase 1の基礎理解（1-2週間）  
**目的**: 既に実装されている検知ロジックを完全に理解する

---

## 📋 目次

1. [粒子データ構造の理解（particle.h）](#1-粒子データ構造の理解particleh)
2. [グリッドセル構造の理解（grid.h）](#2-グリッドセル構造の理解gridh)
3. [検知ロジックの理解（detect_trapped_particles）](#3-検知ロジックの理解detect_trapped_particles)
4. [統合箇所の理解（fix_move_surf.cpp）](#4-統合箇所の理解fix_move_surfcpp)
5. [粒子判定関数の理解（grid_surf.cpp）](#5-粒子判定関数の理解grid_surfcpp)
6. [学習チェックリスト](#学習チェックリスト)

---

## 1. 粒子データ構造の理解（particle.h）

### 1.1 粒子の基本構造

SPARTAでは、粒子は`OnePart`構造体で表現されます。

```cpp
struct SPARTA_ALIGN(16) OnePart {
  int id;                 // 粒子ID（一意の識別子）
  int ispecies;           // 粒子種のインデックス
  int icell;              // 粒子が存在するグリッドセルのインデックス
  int flag;               // マイグレーション状態などに使用
  double x[3];            // 粒子の位置座標 (x, y, z)
  double v[3];            // 粒子の速度ベクトル (vx, vy, vz)
  double erot;            // 回転エネルギー
  double evib;            // 振動エネルギー
  double dtremain;        // 移動タイムステップの残り時間
  double weight;          // 粒子またはセルの重み
};
```

### 1.2 重要なメンバ変数

#### `icell` - セルインデックス
- **意味**: 粒子が現在存在するグリッドセルのインデックス
- **用途**: 粒子がどのセルに属しているかを示す
- **特殊値**: `-1`は削除マーク（粒子が削除対象であることを示す）

#### `x[3]` - 位置座標
- **意味**: 粒子の3次元位置座標
- **用途**: 粒子の位置を表現
- **例**: `x[0] = 1.5, x[1] = 2.3, x[2] = 0.0`（2Dの場合、z=0）

#### `id` - 粒子ID
- **意味**: 粒子の一意の識別子
- **用途**: 粒子を追跡するために使用
- **特徴**: グローバルに一意（全プロセスで一意）

### 1.3 粒子の管理方法

SPARTAでは、粒子はリンクリストで管理されています：

```cpp
Particle::OnePart *particles;  // 粒子配列
int *next;                      // 次の粒子へのリンク
```

**リンクリストの構造**:
- 各セルには`first`（最初の粒子のインデックス）がある
- 各粒子には`next[ip]`（次の粒子のインデックス）がある
- `next[ip] == -1`でリストの終端を示す

**例**:
```
セル0: first = 5
  粒子5: next[5] = 12
  粒子12: next[12] = 23
  粒子23: next[23] = -1  (終端)
```

### 1.4 粒子フラグ（参考）

`particle.h`には粒子フラグの定義は直接含まれていませんが、以下のフラグが使用されます：
- `PKEEP`: 保持する粒子
- `PEXIT`: 出口から出た粒子
- `PENTRY`: 入口から入った粒子
- `PDISCARD`: 破棄する粒子

### 1.5 学習ポイント

✅ **理解すべきこと**:
1. 粒子は`OnePart`構造体で表現される
2. `icell`で粒子がどのセルに属しているかを管理
3. リンクリストでセル内の粒子を管理
4. `icell = -1`は削除マーク

📝 **確認方法**:
```cpp
// 粒子の位置を確認
double *x = particles[ip].x;
printf("Particle %d: pos = (%g, %g, %g)\n", 
       particles[ip].id, x[0], x[1], x[2]);

// 粒子がどのセルに属しているか確認
int icell = particles[ip].icell;
printf("Particle %d is in cell %d\n", particles[ip].id, icell);
```

---

## 2. グリッドセル構造の理解（grid.h）

### 2.1 セルタイプの定義

SPARTAでは、グリッドセルは以下の4つのタイプに分類されます：

```cpp
enum{UNKNOWN, OUTSIDE, INSIDE, OVERLAP};
```

| タイプ | 値 | 意味 | 用途 |
|--------|-----|------|------|
| `UNKNOWN` | 0 | 未分類 | 初期状態やエラー状態 |
| `OUTSIDE` | 1 | 流体領域（壁外） | 粒子が存在できる領域 |
| `INSIDE` | 2 | 固体領域（壁内） | 粒子が存在できない領域 |
| `OVERLAP` | 3 | 境界領域 | 表面要素が含まれるセル |

### 2.2 ChildInfo構造体

各セルには`ChildInfo`構造体で情報が保存されます：

```cpp
struct ChildInfo {
  int count;                // このセル内の粒子数（split cellの場合は0）
  int first;                 // このセル内の最初の粒子のインデックス（-1の場合は粒子なし）
  
  int mask;                  // グリッドグループマスク
  int type;                  // OUTSIDE, INSIDE, OVERLAP, UNKNOWN
  int corner[8];            // コーナーフラグ（2D: 4個、3D: 8個）
  
  double volume;             // セルの体積（流れ体積）
  double weight;             // fnum重み
};
```

### 2.3 セルタイプの判定方法

セルタイプは`grid->set_inout()`関数で分類されます：

1. **OUTSIDE**: セルのすべてのコーナーが壁外
2. **INSIDE**: セルのすべてのコーナーが壁内
3. **OVERLAP**: セルのコーナーが壁内外に混在（表面要素を含む）

### 2.4 セルと粒子の関係

```
セル0 (OUTSIDE):
  count = 5
  first = 10
  粒子10 → 粒子15 → 粒子20 → 粒子25 → 粒子30 → -1

セル1 (INSIDE):
  count = 0
  first = -1
  (粒子なし - 壁内なので粒子は存在しない)

セル2 (OVERLAP):
  count = 2
  first = 40
  粒子40 → 粒子45 → -1
```

### 2.5 学習ポイント

✅ **理解すべきこと**:
1. セルは4つのタイプに分類される
2. `INSIDE`セルは壁内で、通常は粒子が存在しない
3. `OVERLAP`セルは表面要素を含み、粒子が存在する可能性がある
4. `count`と`first`でセル内の粒子を管理

📝 **確認方法**:
```cpp
Grid::ChildInfo *cinfo = grid->cinfo;

// セルのタイプを確認
if (cinfo[icell].type == INSIDE) {
  printf("Cell %d is INSIDE (solid region)\n", icell);
} else if (cinfo[icell].type == OUTSIDE) {
  printf("Cell %d is OUTSIDE (fluid region)\n", icell);
} else if (cinfo[icell].type == OVERLAP) {
  printf("Cell %d is OVERLAP (boundary region)\n", icell);
}

// セル内の粒子数を確認
printf("Cell %d contains %d particles\n", icell, cinfo[icell].count);
```

---

## 3. 検知ロジックの理解（detect_trapped_particles）

### 3.1 関数の目的

`detect_trapped_particles()`は、壁移動後に壁内（固体領域）に取り残された粒子を検知する関数です。

**関数シグネチャ**:
```cpp
void MoveSurf::detect_trapped_particles(int *&trapped_list, int &ntrapped);
```

**引数**:
- `trapped_list`: 検知された粒子のインデックス配列（出力）
- `ntrapped`: 検知された粒子数（出力）

### 3.2 2段階検知システム

この関数は**2段階検知システム**を使用して効率的に粒子を検知します：

#### Stage 1: セルレベルフィルタリング（軽いチェック）

候補セルを絞り込みます：

```cpp
// 候補セルの条件:
// 1. INSIDEタイプのセル
// 2. OVERLAPタイプのセル（表面要素を含む場合）
// 3. 移動した表面を含むセル
```

**コード例** (923-967行目):
```cpp
bool is_candidate = false;

if (cinfo[icell].type == INSIDE) {
  is_candidate = true;
} else if (cinfo[icell].type == OVERLAP) {
  if (cells[icell].nsurf > 0) {
    is_candidate = true;
  }
} else if (cells[icell].nsurf && cells[icell].nsplit >= 1) {
  // 移動した表面を含むセルをチェック
  // pselect配列で移動した点を確認
  ...
}
```

#### Stage 2: 粒子レベル検証（重いチェック）

候補セル内の各粒子を個別に判定します：

```cpp
// 各粒子について outside_surfs() で判定
// outside_surfs() が 0 を返す = 壁内（検知対象）
// outside_surfs() が 1 を返す = 壁外（問題なし）
```

**コード例** (981-1009行目):
```cpp
int ip = cinfo[icell].first;
while (ip >= 0) {
  // 粒子が壁内かどうかを判定
  int is_outside = grid->outside_surfs(icell, particles[ip].x, xcell);
  
  if (is_outside == 0) {
    // 壁内と判定された粒子をリストに追加
    trapped_list[ntrapped_local++] = ip;
  }
  
  // 次の粒子へ
  ip = next[ip];
}
```

### 3.3 処理フロー

```
1. 初期化
   - 粒子が存在しない場合は終了
   - メモリ確保（初期1000個）

2. 全セルをループ
   ├─ Stage 1: 候補セルの判定
   │  ├─ INSIDEセル → 候補
   │  ├─ OVERLAPセル（表面要素あり） → 候補
   │  └─ 移動した表面を含むセル → 候補
   │
   └─ Stage 2: 候補セル内の粒子をチェック
      ├─ セル内の参照点を取得（point_outside_surfs）
      └─ 各粒子について outside_surfs() で判定
         ├─ 壁内（is_outside == 0） → リストに追加
         └─ 壁外（is_outside == 1） → スキップ

3. MPI集約
   - 全プロセスで検知数を集約

4. ログ出力
   - 検知数を出力
   - 最初の10個の粒子の詳細情報を出力
```

### 3.4 重要な関数呼び出し

#### `grid->outside_surfs(icell, x, xcell)`
- **目的**: 粒子が壁外かどうかを判定
- **戻り値**: 1 = 壁外、0 = 壁内
- **引数**:
  - `icell`: セルインデックス
  - `x`: 粒子の位置
  - `xcell`: セル内の参照点（壁外の点）

#### `grid->point_outside_surfs(icell, xcell)`
- **目的**: セル内の参照点（壁外の点）を取得
- **用途**: `outside_surfs()`の引数として使用

### 3.5 メモリ管理

動的メモリ確保を使用：

```cpp
int max_trapped = 1000;
trapped_list = (int *) memory->smalloc(max_trapped*sizeof(int), ...);

// 必要に応じて拡張
if (ntrapped_local >= max_trapped) {
  max_trapped += 1000;
  trapped_list = (int *) memory->srealloc(trapped_list, ...);
}
```

### 3.6 学習ポイント

✅ **理解すべきこと**:
1. 2段階検知システムの動作
2. Stage 1で候補セルを絞り込む理由（効率化）
3. Stage 2で個別粒子を判定する理由（精度確保）
4. `outside_surfs()`の戻り値の意味
5. メモリ管理の方法

📝 **確認方法**:
```cpp
// 検知ロジックの動作を確認
int *trapped_list = NULL;
int ntrapped = 0;
movesurf->detect_trapped_particles(trapped_list, ntrapped);

printf("Detected %d trapped particles\n", ntrapped);
for (int i = 0; i < ntrapped; i++) {
  int ip = trapped_list[i];
  printf("  Particle %d: pos = (%g, %g, %g)\n",
         particles[ip].id,
         particles[ip].x[0], particles[ip].x[1], particles[ip].x[2]);
}
```

---

## 4. 統合箇所の理解（fix_move_surf.cpp）

### 4.1 end_of_step()の処理フロー

`fix_move_surf.cpp`の`end_of_step()`関数は、壁移動処理の全体を統合しています。

**処理順序**:
```
1. 表面移動
   movesurf->move_lines() または movesurf->move_tris()

2. グリッド再構築
   grid->surf2grid()      // 表面要素をグリッドに割り当て
   grid->set_inout()      // セルタイプを分類

3. 検知ロジック
   movesurf->detect_trapped_particles()  ← 追加済み

4. 押し出し処理（Phase 1 月2で実装予定）
   movesurf->push_out_particles()

5. 粒子削除
   movesurf->remove_particles()
```

### 4.2 検知ロジックの呼び出し箇所

**コード例** (261-266行目):
```cpp
// detect trapped particles before removal (for logging and push-out)
int *trapped_list = NULL;
int ntrapped = 0;
if (particle->exist) {
  movesurf->detect_trapped_particles(trapped_list, ntrapped);
}
```

**重要なポイント**:
- `grid->set_inout()`の**後**に呼ばれる
- セルタイプが分類された**後**でないと正しく動作しない
- 削除処理の**前**に呼ばれる（検知してから削除するため）

### 4.3 グリッド再構築の流れ

```cpp
// 1. 表面要素をグリッドに割り当て
grid->surf2grid(1, 0);

// 2. セルタイプを分類
grid->set_inout();

// 3. 検知ロジック実行
movesurf->detect_trapped_particles(trapped_list, ntrapped);
```

### 4.4 学習ポイント

✅ **理解すべきこと**:
1. `end_of_step()`の処理順序
2. 検知ロジックが呼ばれるタイミング
3. グリッド再構築と検知ロジックの関係
4. セルタイプ分類の重要性

📝 **確認方法**:
- `end_of_step()`のコードを読む
- 処理順序を図で整理する
- 実際の実行ログを確認する

---

## 5. 粒子判定関数の理解（grid_surf.cpp）

### 5.1 outside_surfs()関数

**目的**: 粒子が壁外かどうかを判定

**戻り値**:
- `1`: 壁外（流体領域）
- `0`: 壁内（固体領域）

**使用例**:
```cpp
double xcell[3];
grid->point_outside_surfs(icell, xcell);  // 参照点を取得

int is_outside = grid->outside_surfs(icell, particles[ip].x, xcell);
if (is_outside == 0) {
  // 粒子は壁内
}
```

### 5.2 point_outside_surfs()関数

**目的**: セル内の参照点（壁外の点）を取得

**用途**: `outside_surfs()`の引数として使用

**使用例**:
```cpp
double xcell[3];
grid->point_outside_surfs(icell, xcell);
// xcellに壁外の点が設定される
```

### 5.3 EPSSURF定数

```cpp
#define EPSSURF 1.0e-4
```

**意味**: 表面判定のための微小値（安全マージン）

**用途**: 数値誤差を考慮した表面判定

### 5.4 学習ポイント

✅ **理解すべきこと**:
1. `outside_surfs()`の戻り値の意味
2. 参照点（`xcell`）の必要性
3. `EPSSURF`定数の意味

---

## 学習チェックリスト

### ステップ1: 基本構造の理解

- [ ] `particle.h`の`OnePart`構造体を理解した
- [ ] `icell`の意味を理解した
- [ ] リンクリストの構造を理解した
- [ ] `grid.h`の`ChildInfo`構造体を理解した
- [ ] セルタイプ（OUTSIDE, INSIDE, OVERLAP）の意味を理解した

### ステップ2: 検知ロジックの理解

- [ ] `detect_trapped_particles()`の目的を理解した
- [ ] 2段階検知システムの動作を理解した
- [ ] Stage 1（セルレベルフィルタリング）を理解した
- [ ] Stage 2（粒子レベル検証）を理解した
- [ ] `outside_surfs()`の使い方を理解した
- [ ] メモリ管理の方法を理解した

### ステップ3: 統合箇所の理解

- [ ] `fix_move_surf.cpp`の`end_of_step()`の処理フローを理解した
- [ ] 検知ロジックが呼ばれるタイミングを理解した
- [ ] グリッド再構築と検知ロジックの関係を理解した

### ステップ4: 動作確認

- [ ] テストケースを実行した
- [ ] ログ出力を確認した
- [ ] 検知結果を確認した

---

## 実践演習

### 演習1: 粒子データの確認

```cpp
// move_surf.cppのdetect_trapped_particles()内に追加
// 検知された粒子の詳細情報を出力
for (int i = 0; i < ntrapped_local; i++) {
  int ip = trapped_list[i];
  printf("Trapped particle %d:\n", particles[ip].id);
  printf("  Position: (%g, %g, %g)\n",
         particles[ip].x[0], particles[ip].x[1], particles[ip].x[2]);
  printf("  Cell: %d (type: %d)\n", particles[ip].icell, 
         cinfo[particles[ip].icell].type);
}
```

### 演習2: セルタイプの確認

```cpp
// 各セルのタイプを確認
for (int icell = 0; icell < nglocal; icell++) {
  if (cinfo[icell].count > 0) {
    printf("Cell %d: type=%d, count=%d\n", 
           icell, cinfo[icell].type, cinfo[icell].count);
  }
}
```

### 演習3: 2段階検知の動作確認

```cpp
// Stage 1で候補となったセル数をカウント
int n_candidate_cells = 0;
// ... (Stage 1の判定後)
n_candidate_cells++;
printf("Stage 1: %d candidate cells found\n", n_candidate_cells);

// Stage 2で検知された粒子数をカウント
printf("Stage 2: %d trapped particles found\n", ntrapped_local);
```

---

## よくある質問（FAQ）

### Q1: なぜ2段階検知システムを使うのか？

**A**: 計算コストを削減するためです。
- Stage 1（セルレベル）: 軽いチェックで候補を絞り込み
- Stage 2（粒子レベル）: 重いチェックを候補のみに適用

全粒子に対してStage 2を実行すると遅いため、Stage 1で絞り込みます。

### Q2: `outside_surfs()`が0を返すとは？

**A**: 粒子が壁内（固体領域）にあることを意味します。検知対象の粒子です。

### Q3: なぜ参照点（`xcell`）が必要なのか？

**A**: `outside_surfs()`は、粒子から参照点へのレイトレーシングで判定します。参照点は壁外の点である必要があります。

### Q4: `icell = -1`の意味は？

**A**: 粒子が削除対象であることを示します。削除処理で使用されます。

---

## 参考資料

### コード参照
- `src/particle.h`: 粒子データ構造
- `src/grid.h`: グリッドセル構造
- `src/move_surf.cpp`: 検知ロジック実装（895-1035行目）
- `src/fix_move_surf.cpp`: 統合箇所（261-266行目）
- `src/grid_surf.cpp`: 粒子判定関数

### ドキュメント
- `HANDOVER_DOCUMENT.md`: 包括的な引き継ぎ資料
- `DETECTION_IMPLEMENTATION.md`: 検知ロジック実装報告
- `DETECTION_VERIFICATION.md`: 検証結果の詳細

---

**最終更新**: 2024年12月26日  
**ステータス**: Phase 1基礎理解用勉強資料完成
