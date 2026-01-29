# 押し出し処理実装 引き継ぎ資料

## 📋 現在の状態

### ✅ 完了している作業
1. **検知ロジックの実装**: `detect_trapped_particles()` が実装済み
2. **削除ロジックの修正**: `remove_particles()` が個別粒子判定に対応
3. **検知数と削除数の整合性**: ほぼ一致（高密度時は稀に不一致）
4. **OVERLAPセルの対応**: 検知・削除の両方で対応済み
5. **Split cellの対応**: 削除マークの保持・再適用に対応

### ❌ 未実装の機能
- **押し出し処理**: `push_out_particles()` は未実装
- **表面速度計算**: `compute_surf_velocity()` は未実装
- **幾何計算関数**: 点から表面要素への最短距離計算（参考関数は既存）

## 🎯 実装すべき機能

### 1. 押し出し処理の全体フロー

**関数**: `MoveSurf::push_out_particles()`

**処理フロー**:
1. `detect_trapped_particles()` で検知された粒子のリストを使用
2. 各粒子について押し出し処理を実行
3. 押し出し成功した粒子は削除しない（`icell` を維持または更新）
4. 押し出し失敗した粒子は削除する（`icell = -1` を設定）

**実装場所**: `src/move_surf.cpp`

### 2. 5ステップアルゴリズム

#### Step 1: Cell Selection
**目的**: 粒子が存在するセルの表面要素を取得

**実装**:
```cpp
// 粒子の現在のセルを取得
int icell = particles[ip].icell;
int nsurf = cells[icell].nsurf;
surfint *csurfs = cells[icell].csurfs;

// このセルに含まれる表面要素を取得
for (int m = 0; m < nsurf; m++) {
  int isurf = csurfs[m];
  // この表面要素を候補として検討
}
```

#### Step 2: Velocity Filter
**目的**: 押し出し方向の表面要素のみを候補とする

**実装**:
```cpp
// 表面速度を計算（並進・回転を考慮）
double vwall[3];
compute_surf_velocity(isurf, particles[ip].x, vwall);

// 法線ベクトルを取得
double *norm;
if (dim == 3) {
  norm = surf->tris[isurf].norm;
} else {
  norm = surf->lines[isurf].norm;
}

// 法線方向の速度成分を計算
double V_n = vwall[0]*norm[0] + vwall[1]*norm[1] + vwall[2]*norm[2];

// 押し出し方向の判定（V_n > 0）
if (V_n <= 0.0) continue;  // この表面要素はスキップ
```

#### Step 3: Projection Check
**目的**: 射影パラメータ $0 \le t \le 1$ で判定

**実装**:
```cpp
// 点から線分/三角形への最短距離と射影パラメータを計算
double dist, t;
double proj_pos[3];

if (dim == 2) {
  // 2D: 点から線分への最短距離
  compute_distance_to_line(particles[ip].x, isurf, dist, t, proj_pos);
} else {
  // 3D: 点から三角形への最短距離
  compute_distance_to_triangle(particles[ip].x, isurf, dist, t, proj_pos);
}

// 射影パラメータの範囲チェック
// 注意: ここでの t は射影パラメータ（0-1）とは異なる可能性がある
// geometry.cppの関数を参考に実装する必要がある
```

**注意**: `geometry.cpp` の `distsq_point_line()` と `distsq_point_tri()` を参考にするが、これらの関数は射影パラメータ `t` を返さない。新しく実装する必要がある。

#### Step 4: Re-positioning
**目的**: 最短距離の面に安全マージンを加えて再配置

**実装**:
```cpp
// 最短距離の表面を選択（複数候補がある場合は最短距離のものを選択）
// 安全マージンの計算
double epsilon = max(EPSSURF * min_edge, grid->cell_epsilon);
// EPSSURF = 1.0e-4 (grid_surf.cpp で定義)

// 新しい位置を計算
double new_pos[3];
new_pos[0] = proj_pos[0] + epsilon * norm[0];
new_pos[1] = proj_pos[1] + epsilon * norm[1];
new_pos[2] = proj_pos[2] + epsilon * norm[2];

// 粒子位置を更新
particles[ip].x[0] = new_pos[0];
particles[ip].x[1] = new_pos[1];
particles[ip].x[2] = new_pos[2];
```

#### Step 5: Velocity Update
**目的**: CLLモデル（または他の表面衝突モデル）を使用して速度を更新

**実装**:
```cpp
// 表面要素から衝突モデルを取得
Surf::Tri *tri = NULL;
Surf::Line *line = NULL;
int isc = -1;
double *norm = NULL;

if (dim == 3) {
  tri = &surf->tris[isurf];
  norm = tri->norm;
  isc = tri->isc;
} else {
  line = &surf->lines[isurf];
  norm = line->norm;
  isc = line->isc;
}

// CLLモデル（または他の表面衝突モデル）を使用して速度を更新
// 注意: 通常の衝突とは異なり、粒子は既に表面外に配置されている
// したがって、dtremain は実際には使用されない（ダミー値でOK）
double dtremain = update->dt;  // ダミー値（実際には使用されない）
int reaction = 0;
Particle::OnePart *jp = surf->sc[isc]->collide(
    &particles[ip], dtremain, isurf, norm, -1, reaction);
// 第4引数: isurf (表面要素インデックス)
// 第5引数: -1 (表面反応モデルなし)
// 第6引数: reaction (出力、反応インデックス)

// jp は新しい粒子（化学反応で生成された場合）
// 通常は NULL を返す（無視してOK）
```

**参考コード**: `update.cpp:918-923`

### 3. 必要な関数の実装

#### 関数1: `compute_surf_velocity()`
**目的**: 表面要素の速度を計算（並進・回転を考慮）

**シグネチャ**:
```cpp
void MoveSurf::compute_surf_velocity(int isurf, double *x, double *vwall);
```

**実装方針**:
- **並進の場合**: `vwall = (vx, vy, vz)` （固定値）
  - `fix_move_surf` から速度パラメータを取得する必要がある
  - または、`surf_collide` から速度パラメータを取得
- **回転の場合**: `vwall = ω × (x - p)` （位置依存）
  - `vwall[0] = wy*(x[2]-pz) - wz*(x[1]-py);`
  - `vwall[1] = wz*(x[0]-px) - wx*(x[2]-pz);`
  - `vwall[2] = wx*(x[1]-py) - wy*(x[0]-px);`

**注意点**:
- `move_surf` は現在、表面の移動のみを処理（変位ベース）
- 速度情報は `fix_move_surf` に保持されている可能性がある
- または、`surf_collide` から速度パラメータを取得する必要がある
- **重要**: `surf_collide_cll` には速度パラメータ（vx, vy, vz, wx, wy, wz）があるが、これは `surf_collide` オブジェクトのパラメータであり、表面要素ごとではない

**解決策**:
- 表面要素に割り当てられた `surf_collide` オブジェクトから速度パラメータを取得
- または、`fix_move_surf` から速度パラメータを取得（現在の実装では困難）

#### 関数2: `compute_distance_to_line()` (2D)
**目的**: 点から線分への最短距離と射影パラメータを計算

**シグネチャ**:
```cpp
void MoveSurf::compute_distance_to_line(
    double *x, int isurf, double &dist, double &t, double *proj_pos);
```

**実装方針**:
- 線分の端点: `surf->lines[isurf].p1`, `surf->lines[isurf].p2`
- `geometry.cpp` の `distsq_point_line()` を参考にする
- 射影パラメータ `t` を計算: `t = dot(a, b) / lensq(b)`
- `t` を `[0,1]` にクランプ
- 最短距離と射影位置を計算

**参考コード**: `geometry.cpp:1251-1282`

#### 関数3: `compute_distance_to_triangle()` (3D)
**目的**: 点から三角形への最短距離と射影パラメータを計算

**シグネチャ**:
```cpp
void MoveSurf::compute_distance_to_triangle(
    double *x, int isurf, double &dist, double &t, double *proj_pos);
```

**実装方針**:
- 三角形の頂点: `surf->tris[isurf].p1`, `p2`, `p3`
- 法線ベクトル: `surf->tris[isurf].norm`
- `geometry.cpp` の `distsq_point_tri()` を参考にする
- 点から三角形平面への射影を計算
- 射影点が三角形内にあるか判定（重心座標を使用）
- 三角形外の場合は、エッジまたは頂点への最短距離を計算
- **注意**: 射影パラメータ `t` の定義が曖昧（重心座標を使用する可能性がある）

**参考コード**: `geometry.cpp:1289-1344`

#### 関数4: `push_out_particles()`
**目的**: 検知された粒子を押し出す

**シグネチャ**:
```cpp
int MoveSurf::push_out_particles(int *trapped_list, int ntrapped);
```

**戻り値**: 押し出しに成功した粒子数

**実装方針**:
- `trapped_list` に含まれる各粒子について押し出し処理を実行
- 5ステップアルゴリズムを適用
- 押し出し成功した粒子は削除マークを付けない
- 押し出し失敗した粒子は削除マークを付ける（`icell = -1`）
- 押し出し後に `icell` を更新する必要がある

### 4. fix_move_surf.cpp での統合

**現在の処理フロー** (`fix_move_surf.cpp:180-208`):
```cpp
// detect trapped particles before removal (for logging)
int *trapped_list = NULL;
int ntrapped = 0;
if (particle->exist) {
  movesurf->detect_trapped_particles(trapped_list, ntrapped);
}

// remove particles as needed due to surface move
if (particle->exist) ndeleted = movesurf->remove_particles();

// Log deletion count for comparison with detection
if (me == 0) {
  bigint ntrapped_global;
  bigint ntrapped_big = ntrapped;
  MPI_Allreduce(&ntrapped_big, &ntrapped_global, 1, 
                MPI_SPARTA_BIGINT, MPI_SUM, world);
  if (screen) {
    fprintf(screen, "Deleted " BIGINT_FORMAT " particles (detected: " BIGINT_FORMAT ")\n",
            ndeleted, ntrapped_global);
  }
  if (logfile) {
    fprintf(logfile, "Deleted " BIGINT_FORMAT " particles (detected: " BIGINT_FORMAT ")\n",
            ndeleted, ntrapped_global);
  }
}

// free trapped_list if allocated
if (trapped_list) memory->sfree(trapped_list);
```

**変更後の処理フロー**:
```cpp
// detect trapped particles
int *trapped_list = NULL;
int ntrapped = 0;
if (particle->exist) {
  movesurf->detect_trapped_particles(trapped_list, ntrapped);
}

// push out trapped particles
int npushed = 0;
if (particle->exist && ntrapped > 0) {
  npushed = movesurf->push_out_particles(trapped_list, ntrapped);
}

// remove particles that could not be pushed out
if (particle->exist) ndeleted = movesurf->remove_particles();

// Log results
if (me == 0) {
  bigint ntrapped_global, npushed_global;
  bigint ntrapped_big = ntrapped;
  bigint npushed_big = npushed;
  MPI_Allreduce(&ntrapped_big, &ntrapped_global, 1, 
                MPI_SPARTA_BIGINT, MPI_SUM, world);
  MPI_Allreduce(&npushed_big, &npushed_global, 1, 
                MPI_SPARTA_BIGINT, MPI_SUM, world);
  if (screen) {
    fprintf(screen, "Pushed out " BIGINT_FORMAT " particles, deleted " BIGINT_FORMAT 
            " particles (detected: " BIGINT_FORMAT ")\n",
            npushed_global, ndeleted, ntrapped_global);
  }
  if (logfile) {
    fprintf(logfile, "Pushed out " BIGINT_FORMAT " particles, deleted " BIGINT_FORMAT 
            " particles (detected: " BIGINT_FORMAT ")\n",
            npushed_global, ndeleted, ntrapped_global);
  }
}

// free trapped_list if allocated
if (trapped_list) memory->sfree(trapped_list);
```

## ⚠️ 重要な注意点

### 1. セル再割り当て（icell の更新）
**問題**: 押し出し後に粒子位置を変更すると、`icell` が不整合になる

**解決策**:
```cpp
// 押し出し後に icell を更新
// 参考: grid_adapt.cpp:154-169

cellint parentID = cells[icell].id;
int plevel = cells[icell].level;
double *plo = cells[icell].lo;
double *phi = cells[icell].hi;

int new_icell = grid->id_find_child(parentID, plevel, plo, phi, particles[ip].x);

// split cell の場合は sub cell を特定
if (cells[new_icell].nsplit > 1) {
  if (dim == 3) {
    new_icell = update->split3d(new_icell, particles[ip].x);
  } else {
    new_icell = update->split2d(new_icell, particles[ip].x);
  }
}

particles[ip].icell = new_icell;
```

**参考コード**: `grid_adapt.cpp:146-185`

### 2. リンクリストの更新
**問題**: `icell` が変更されると、リンクリストが不整合になる

**解決策**:
```cpp
// 押し出し処理中は、リンクリストを手動で更新する必要がある
// または、押し出し処理後に particle->sort() を呼ぶ

// 推奨: 押し出し処理後に particle->sort() を呼ぶ
// これにより、すべての粒子が正しいセルに割り当てられる
particle->sort();
```

### 3. 速度更新の必要性
**理由**: 壁が粒子より速い場合、速度を更新しないと再取り込まれる

**実装**: Step 5でCLLモデルを使用

**注意**: 
- CLLモデルは移動壁を考慮している（`surf_collide_cll.cpp` を参照）
- `trflag` が設定されている場合、並進または回転速度が自動的に追加される

### 4. 高密度時の処理時間
**問題**: 粒子数が多い場合、処理時間が増加する可能性

**対策**: 
- 効率的なアルゴリズムの実装
- 必要に応じて最適化
- 並列処理の考慮（現在は各プロセスで独立に処理）

### 5. 表面速度の取得方法
**問題**: `move_surf` には速度情報がない（変位ベース）

**解決策のオプション**:
1. **`surf_collide` から速度パラメータを取得**
   - `surf_collide_cll` には `vx, vy, vz, wx, wy, wz` がある
   - しかし、これらは表面要素ごとではなく、`surf_collide` オブジェクトごと
   - 表面要素に割り当てられた `surf_collide` オブジェクトを取得する必要がある

2. **`fix_move_surf` から速度パラメータを取得**
   - `fix_move_surf` は `move_surf` のインスタンスを持つ
   - しかし、`move_surf` は変位ベースで速度情報を保持していない
   - `fix_move_surf` に速度情報を追加する必要がある可能性

3. **時間ステップから速度を計算**
   - `fix_move_surf` は `nevery` ごとに移動を実行
   - 変位から速度を計算: `v = delta / (nevery * dt)`
   - しかし、回転の場合、速度は位置依存なので、この方法は不適切

**推奨**: オプション1（`surf_collide` から速度パラメータを取得）を試す

## 📁 関連ファイル

### 実装するファイル
- `src/move_surf.h`: 関数宣言を追加
- `src/move_surf.cpp`: 関数実装
- `src/fix_move_surf.cpp`: 呼び出し箇所

### 参考になるファイル
- `src/surf_collide_cll.cpp`: CLLモデルの実装（速度更新の参考）
- `src/surf_collide_cll.h`: CLLモデルの定義（速度パラメータの参照）
- `src/update.cpp`: 通常の衝突処理（参考）
- `src/geometry.cpp`: 距離計算関数（参考）
- `src/grid_adapt.cpp`: セル再割り当て（参考）
- `src/grid_surf.cpp`: EPSSURF の定義（参考）

## 📚 参考になるコード

### 通常の衝突処理（update.cpp:918-923）
```cpp
if (DIM == 3)
  jpart = surf->sc[tri->isc]->
    collide(ipart,dtremain,minsurf,tri->norm,tri->isr,reaction);
if (DIM != 3)
  jpart = surf->sc[line->isc]->
    collide(ipart,dtremain,minsurf,line->norm,line->isr,reaction);
```

### セル再割り当て（grid_adapt.cpp:154-169）
```cpp
ichild = id_find_child(parentID,plevel,plo,phi,p->x);

if (ichild < 0) {
  printf("BAD CHILD %d: %d " CELLINT_FORMAT ": %g %g %g\n",
         me,icell,parentID,p->x[0],p->x[1],p->x[2]);
  error->one(FLERR,"Adapt refine particle could not be mapped to child cell");
}

if (cells[ichild].nsplit > 1) {
  if (dim == 3) ichild = update->split3d(ichild,p->x);
  else ichild = update->split2d(ichild,p->x);
}
```

### 距離計算（geometry.cpp:1251-1282）
```cpp
double distsq_point_line(double *x, double *p1, double *p2)
{
  double a[3],b[3],c[3];
  MathExtra::sub3(x,p1,a);
  MathExtra::sub3(p2,p1,b);
  
  double alpha = MathExtra::dot3(a,b)/MathExtra::lensq3(b);
  
  if (alpha >= 1.0) MathExtra::sub3(x,p2,c);
  else if (alpha > 0.0) {
    a[0] = p1[0] + alpha*b[0];
    a[1] = p1[1] + alpha*b[1];
    a[2] = p1[2] + alpha*b[2];
    MathExtra::sub3(x,a,c);
  } else MathExtra::sub3(x,p1,c);
  
  return MathExtra::lensq3(c);
}
```

### CLLモデルの速度更新（surf_collide_cll.cpp）
```cpp
// 並進の場合
if (tflag) {
  vxdelta = vx; vydelta = vy; vzdelta = vz;
  // ... (normal velocity component correction) ...
}

// 回転の場合
else if (rflag) {
  double *x = p->x;
  vxdelta = wy*(x[2]-pz) - wz*(x[1]-py);
  vydelta = wz*(x[0]-px) - wx*(x[2]-pz);
  vzdelta = wx*(x[1]-py) - wy*(x[0]-px);
  double dot = vxdelta*norm[0] + vydelta*norm[1] + vzdelta*norm[2];
  vxdelta -= dot*norm[0];
  vydelta -= dot*norm[1];
  vzdelta -= dot*norm[2];
}

v[0] = vperp*norm[0] + vtan1*tangent1[0] + vtan2*tangent2[0] + vxdelta;
v[1] = vperp*norm[1] + vtan1*tangent1[1] + vtan2*tangent2[1] + vydelta;
v[2] = vperp*norm[2] + vtan1*tangent1[2] + vtan2*tangent2[2] + vzdelta;
```

## ✅ 実装チェックリスト

### Phase 1: 基本的な押し出し処理（位置のみ更新）
- [ ] `move_surf.h` に関数宣言を追加
  - [ ] `void compute_surf_velocity(int isurf, double *x, double *vwall);`
  - [ ] `void compute_distance_to_line(double *x, int isurf, double &dist, double &t, double *proj_pos);`
  - [ ] `void compute_distance_to_triangle(double *x, int isurf, double &dist, double &t, double *proj_pos);`
  - [ ] `int push_out_particles(int *trapped_list, int ntrapped);`
- [ ] `move_surf.cpp` に関数を実装
  - [ ] Step 1: Cell Selection
  - [ ] Step 2: Velocity Filter
  - [ ] Step 3: Projection Check
  - [ ] Step 4: Re-positioning
  - [ ] セル再割り当ての実装
- [ ] `fix_move_surf.cpp` で呼び出しを統合
- [ ] リンクリストの更新（`particle->sort()`）
- [ ] エラーチェック（`fix grid/check`）
- [ ] テスト（低密度で動作確認）

### Phase 2: 速度更新を追加
- [ ] Step 5: Velocity Update を実装
  - [ ] CLLモデルの使用
  - [ ] 表面要素から `surf_collide` オブジェクトを取得
  - [ ] `collide()` 関数の呼び出し
- [ ] テスト（低密度で動作確認）

### Phase 3: 完全な実装
- [ ] 高密度でのテスト
- [ ] エラーチェック
- [ ] パフォーマンス最適化
- [ ] ドキュメント更新

## 🔧 実装順序の推奨

1. **Phase 1**: 基本的な押し出し処理（位置のみ更新、速度は変更しない）
   - Step 1-4 を実装
   - セル再割り当てを実装
   - テスト（低密度）

2. **Phase 2**: 速度更新を追加
   - Step 5 を実装（CLLモデルの使用）
   - テスト（低密度）

3. **Phase 3**: 完全な実装
   - 高密度でのテスト
   - エラーチェック
   - パフォーマンス最適化

## 🐛 トラブルシューティング

### 問題1: icell が正しく更新されない
- **症状**: `fix grid/check` でエラー
- **解決策**: `grid->id_find_child()` と `update->split2d/split3d()` を使用
- **参考**: `grid_adapt.cpp:154-169`

### 問題2: リンクリストが不整合になる
- **症状**: 粒子が見つからない、または重複
- **解決策**: `particle->sort()` を呼ぶか、手動でリンクリストを更新
- **参考**: `particle->sort()` の実装

### 問題3: 速度を更新しても再取り込まれる
- **症状**: 押し出し後に粒子が再び壁内に戻る
- **解決策**: CLLモデルの移動壁パラメータが正しく設定されているか確認
- **参考**: `surf_collide_cll.cpp` の `trflag` と速度パラメータ

### 問題4: 表面速度が取得できない
- **症状**: `compute_surf_velocity()` で速度がゼロまたは不正
- **解決策**: `surf_collide` から速度パラメータを取得する方法を確認
- **参考**: `surf_collide_cll.h` の `vx, vy, vz, wx, wy, wz`

### 問題5: 射影パラメータの計算が難しい
- **症状**: Step 3 で適切な表面要素を選択できない
- **解決策**: `geometry.cpp` の関数を参考に、射影パラメータを正しく計算
- **参考**: `geometry.cpp:1251-1344`

## 📖 関連ドキュメント

- `HANDOVER_DOCUMENT.md`: プロジェクト全体の計画
- `IMPLEMENTATION_PHASE1.md`: Phase 1 の実装計画
- `PLAN_REVIEW.md`: 計画のレビュー

## 🔗 重要な定数と定義

- `EPSSURF`: 1.0e-4 (`grid_surf.cpp:44`)
- `UNKNOWN, OUTSIDE, INSIDE, OVERLAP`: セルタイプ（`grid.cpp` など）
- `PKEEP, PEXIT, PENTRY, PDISCARD`: 粒子フラグ（`particle.h`）

## 💡 実装のヒント

1. **デバッグログの追加**: 押し出し処理の各ステップでログを出力して、問題を特定しやすくする
2. **段階的な実装**: Step 1 から順番に実装し、各ステップでテストする
3. **既存コードの再利用**: `geometry.cpp` の関数を参考にして、新しい関数を実装する
4. **エラーハンドリング**: 不正な状態（例: 粒子が見つからない、セルが存在しない）を適切に処理する

---

**最終更新**: 2024年（実装開始時点）  
**次のステップ**: Phase 1 の実装を開始する
