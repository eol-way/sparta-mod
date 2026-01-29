# SPARTA改造プロジェクト：関連ファイル一覧

**作成日**: 2024年12月26日  
**目的**: 改造計画に関連するh, cppファイルをまとめたドキュメント

---

## 📋 目次

1. [実装・修正ファイル（直接変更）](#実装修正ファイル直接変更)
2. [参考ファイル（参照・参考にする）](#参考ファイル参照参考にする)
3. [ファイル別の役割と実装箇所](#ファイル別の役割と実装箇所)
4. [実装フェーズ別ファイル分類](#実装フェーズ別ファイル分類)

---

## 実装・修正ファイル（直接変更）

### 1. `src/move_surf.h`
**役割**: MoveSurfクラスの定義と関数宣言  
**変更内容**:
- `detect_trapped_particles(int *&, int &)` 関数宣言を追加（Phase 1 月1完了）
- `push_out_particles(int *trapped_list, int ntrapped)` 関数宣言を追加予定（Phase 1 月2）
- `compute_surf_velocity(int isurf, double *x, double *vwall)` 関数宣言を追加予定（Phase 1 月2）
- `compute_distance_to_line(double *x, int isurf, double &dist, double &t, double *proj_pos)` 関数宣言を追加予定（Phase 1 月2）
- `compute_distance_to_triangle(double *x, int isurf, double &dist, double &t, double *proj_pos)` 関数宣言を追加予定（Phase 1 月2）

**参照ドキュメント**: 
- `HANDOVER_DOCUMENT.md` (161行目)
- `DETECTION_IMPLEMENTATION.md` (7行目)
- `SESSION_SUMMARY.md` (94行目)

---

### 2. `src/move_surf.cpp`
**役割**: MoveSurfクラスの実装  
**変更内容**:
- `detect_trapped_particles()` 関数実装（Phase 1 月1完了、約895行目以降）
- `push_out_particles()` 関数実装予定（Phase 1 月2）
- `compute_surf_velocity()` 関数実装予定（Phase 1 月2）
- `compute_distance_to_line()` 関数実装予定（Phase 1 月2）
- `compute_distance_to_triangle()` 関数実装予定（Phase 1 月2）

**主要機能**:
- 2段階検知システム（セルレベル + 粒子レベル）
- 押し出しアルゴリズム（5ステップ）
- 表面速度計算（回転・並進対応）

**参照ドキュメント**: 
- `HANDOVER_DOCUMENT.md` (167行目, 486行目)
- `DETECTION_IMPLEMENTATION.md` (10行目)
- `SESSION_SUMMARY.md` (95行目)
- `PUSH_OUT_HANDOVER.md` (29行目, 405行目)

---

### 3. `src/fix_move_surf.cpp`
**役割**: Fixクラスのラッパー、壁移動処理の統合  
**変更内容**:
- `end_of_step()` 内で検知ロジックを呼び出し（Phase 1 月1完了、180-194行目）
- 押し出し処理の統合予定（Phase 1 月2）
- 削除処理を押し出し処理に置き換え予定（Phase 1 月2）

**処理フロー**:
1. 表面移動 (`movesurf->move_lines/tris()`)
2. グリッド再構築 (`grid->surf2grid()`, `grid->set_inout()`)
3. 検知ロジック (`movesurf->detect_trapped_particles()`) ← 追加済み
4. 押し出し処理 (`movesurf->push_out_particles()`) ← 月2で追加予定
5. 粒子削除 (`movesurf->remove_particles()`) ← 月2で置き換え予定

**参照ドキュメント**: 
- `HANDOVER_DOCUMENT.md` (185行目, 487行目)
- `DETECTION_IMPLEMENTATION.md` (35行目)
- `SESSION_SUMMARY.md` (246行目)
- `PUSH_OUT_HANDOVER.md` (244行目, 406行目)

---

### 4. `src/fix_move_surf.h`
**役割**: FixMoveSurfクラスの定義  
**変更内容**: 変更なし（言及のみ）

**参照ドキュメント**: 
- `HANDOVER_DOCUMENT.md` (210行目)

---

## 参考ファイル（参照・参考にする）

### 5. `src/surf_collide_cll.cpp`
**役割**: CLL（Cercignani-Lampis-Lord）表面衝突モデルの実装  
**参考箇所**:
- 速度更新の実装（押し出し後の速度更新に使用）
- 移動壁を考慮した速度計算（`trflag`による並進・回転速度の追加）
- `collide()` 関数の呼び出し方法

**参照ドキュメント**: 
- `PUSH_OUT_HANDOVER.md` (409行目, 466行目, 554行目)

---

### 6. `src/surf_collide_cll.h`
**役割**: CLLモデルの定義と速度パラメータ  
**参考箇所**:
- 速度パラメータ（`vx, vy, vz, wx, wy, wz`）の定義
- 並進・回転フラグ（`tflag, rflag`）の定義
- 回転中心（`px, py, pz`）の定義

**参照ドキュメント**: 
- `PUSH_OUT_HANDOVER.md` (410行目, 559行目)

---

### 7. `src/surf_collide_diffuse.h`
**役割**: 拡散反射モデルの定義  
**参考箇所**:
- 移動壁の速度計算機能（既に実装済み）
- 速度パラメータ（`vx, vy, vz, wx, wy, wz`）の参照

**参照ドキュメント**: 
- `PLAN_REVIEW.md` (37行目)

---

### 8. `src/update.cpp`
**役割**: 通常の衝突処理  
**参考箇所**:
- 918-923行目: 通常の衝突処理の実装
- `surf->sc[isc]->collide()` の呼び出し方法
- 衝突モデルの使用方法

**参照ドキュメント**: 
- `PUSH_OUT_HANDOVER.md` (156行目, 411行目, 418行目)

---

### 9. `src/geometry.cpp`
**役割**: 幾何計算関数（距離計算など）  
**参考箇所**:
- 1251-1282行目: `distsq_point_line()` - 点から線分への距離計算（2D用）
- 1289-1344行目: `distsq_point_tri()` - 点から三角形への距離計算（3D用）
- 射影パラメータの計算方法（ただし、既存関数は射影パラメータを返さないため、新規実装が必要）

**参照ドキュメント**: 
- `PUSH_OUT_HANDOVER.md` (93行目, 96行目, 198行目, 203行目, 217行目, 223行目, 412行目, 444行目, 563行目, 564行目)

---

### 10. `src/grid_adapt.cpp`
**役割**: グリッド適応処理（セル再割り当て）  
**参考箇所**:
- 146-185行目: セル再割り当ての実装
- 154-169行目: `id_find_child()` と `split2d/split3d()` の使用例
- 押し出し後のセル割り当て更新に使用

**参照ドキュメント**: 
- `PUSH_OUT_HANDOVER.md` (327行目, 348行目, 413行目, 428行目, 544行目)

---

### 11. `src/grid_surf.cpp`
**役割**: グリッドと表面の相互作用、粒子判定関数  
**参考箇所**:
- 44行目: `EPSSURF` 定数の定義（1.0e-4）
- `outside_surfs()` 関数: 個別粒子の壁内外判定
- `point_outside_surfs()` 関数: セル内の参照点取得

**参照ドキュメント**: 
- `HANDOVER_DOCUMENT.md` (488行目)
- `PLAN_REVIEW.md` (114行目)
- `PUSH_OUT_HANDOVER.md` (414行目, 574行目)

---

### 12. `src/grid.h`
**役割**: Gridクラスの定義、セルタイプ定義  
**参考箇所**:
- 150-165行目: `ChildInfo` 構造体（セルタイプ定義）
- セルタイプ: `UNKNOWN`, `OUTSIDE`, `INSIDE`, `OVERLAP`

**参照ドキュメント**: 
- `SESSION_SUMMARY.md` (247行目)
- `PUSH_OUT_HANDOVER.md` (575行目)

---

### 13. `src/grid.cpp`
**役割**: Gridクラスの実装  
**参考箇所**:
- 1428行目以降: `set_inout()` 関数（セルタイプの分類）
- セルタイプ判定のロジック

**参照ドキュメント**: 
- `SESSION_SUMMARY.md` (248行目)
- `PUSH_OUT_HANDOVER.md` (575行目)

---

### 14. `src/particle.h`
**役割**: 粒子データ構造の定義  
**参考箇所**:
- 粒子フラグ: `PKEEP`, `PEXIT`, `PENTRY`, `PDISCARD`
- 粒子データ構造の理解

**参照ドキュメント**: 
- `HANDOVER_DOCUMENT.md` (198行目)
- `PUSH_OUT_HANDOVER.md` (576行目)

---

## ファイル別の役割と実装箇所

### Phase 1 月1（検知ロジック）✅ 完了

| ファイル | 役割 | 実装箇所 | 状態 |
|---------|------|---------|------|
| `move_surf.h` | 関数宣言 | `detect_trapped_particles()` | ✅ 完了 |
| `move_surf.cpp` | 関数実装 | 約895行目以降 | ✅ 完了 |
| `fix_move_surf.cpp` | 統合 | 180-194行目 | ✅ 完了 |

---

### Phase 1 月2（押し出しロジック）⏭️ 未着手

| ファイル | 役割 | 実装箇所 | 状態 |
|---------|------|---------|------|
| `move_surf.h` | 関数宣言追加 | `push_out_particles()`, `compute_surf_velocity()`, `compute_distance_to_line()`, `compute_distance_to_triangle()` | ⏭️ 未着手 |
| `move_surf.cpp` | 関数実装 | 押し出しアルゴリズム（5ステップ） | ⏭️ 未着手 |
| `fix_move_surf.cpp` | 統合 | 押し出し処理の呼び出し、削除処理の置き換え | ⏭️ 未着手 |

**参考ファイル**:
- `surf_collide_cll.cpp/h`: 速度更新の参考
- `geometry.cpp`: 距離計算の参考
- `grid_adapt.cpp`: セル再割り当ての参考
- `update.cpp`: 衝突処理の参考

---

### Phase 1 月3（並列化対応）⏭️ 未着手

| ファイル | 役割 | 実装箇所 | 状態 |
|---------|------|---------|------|
| `fix_move_surf.cpp` | MPI通信処理 | 粒子のプロセス間移動 | ⏭️ 未着手 |
| `move_surf.cpp` | 並列化対応 | Ghost cell更新 | ⏭️ 未着手 |

---

## 実装フェーズ別ファイル分類

### フェーズ1: 検知ロジック（✅ 完了）

**直接変更ファイル**:
1. `src/move_surf.h` - 関数宣言追加
2. `src/move_surf.cpp` - 検知ロジック実装
3. `src/fix_move_surf.cpp` - 統合

**参考ファイル**:
- `src/grid_surf.cpp` - 粒子判定関数
- `src/grid.h/cpp` - セルタイプ定義・設定

---

### フェーズ2: 押し出しロジック（⏭️ 未着手）

**直接変更ファイル**:
1. `src/move_surf.h` - 関数宣言追加
2. `src/move_surf.cpp` - 押し出しアルゴリズム実装
3. `src/fix_move_surf.cpp` - 統合・置き換え

**参考ファイル**:
- `src/surf_collide_cll.cpp/h` - 速度更新の参考
- `src/geometry.cpp` - 距離計算の参考
- `src/grid_adapt.cpp` - セル再割り当ての参考
- `src/update.cpp` - 衝突処理の参考
- `src/grid_surf.cpp` - EPSSURF定義

---

### フェーズ3: 並列化対応（⏭️ 未着手）

**直接変更ファイル**:
1. `src/fix_move_surf.cpp` - MPI通信処理
2. `src/move_surf.cpp` - 並列化対応

**参考ファイル**:
- 既存のMPI通信処理を参考にする

---

## 重要な定数と定義

### 定数
- `EPSSURF`: 1.0e-4 (`grid_surf.cpp:44`)
- セルタイプ: `UNKNOWN`, `OUTSIDE`, `INSIDE`, `OVERLAP` (`grid.h`)
- 粒子フラグ: `PKEEP`, `PEXIT`, `PENTRY`, `PDISCARD` (`particle.h`)

---

## 実装チェックリスト

### Phase 1 月1（検知ロジック）✅ 完了
- [x] `move_surf.h` に関数宣言を追加
- [x] `move_surf.cpp` に関数を実装
- [x] `fix_move_surf.cpp` で呼び出しを統合
- [x] テスト（低密度で動作確認）

### Phase 1 月2（押し出しロジック）⏭️ 未着手
- [ ] `move_surf.h` に関数宣言を追加
  - [ ] `push_out_particles()`
  - [ ] `compute_surf_velocity()`
  - [ ] `compute_distance_to_line()`
  - [ ] `compute_distance_to_triangle()`
- [ ] `move_surf.cpp` に関数を実装
  - [ ] Step 1: Cell Selection
  - [ ] Step 2: Velocity Filter
  - [ ] Step 3: Projection Check
  - [ ] Step 4: Re-positioning
  - [ ] Step 5: Velocity Update
  - [ ] セル再割り当ての実装
- [ ] `fix_move_surf.cpp` で呼び出しを統合
- [ ] リンクリストの更新（`particle->sort()`）
- [ ] テスト（低密度で動作確認）

### Phase 1 月3（並列化対応）⏭️ 未着手
- [ ] 粒子のプロセス間移動処理
- [ ] Ghost cell更新の統合
- [ ] メモリリークチェック
- [ ] ベンチマークテスト（BM-1）

---

## 参考ドキュメント

詳細な実装計画や技術的説明については、以下のドキュメントを参照してください：

1. **`HANDOVER_DOCUMENT.md`**: 包括的な引き継ぎ資料
2. **`PUSH_OUT_HANDOVER.md`**: 押し出し処理の詳細な実装計画
3. **`PLAN_REVIEW.md`**: 計画レビューと改善提案
4. **`IMPLEMENTATION_PHASE1.md`**: Phase 1の実装詳細計画
5. **`DETECTION_IMPLEMENTATION.md`**: 検知ロジック実装報告
6. **`DETECTION_VERIFICATION.md`**: 検証結果の詳細
7. **`SESSION_SUMMARY.md`**: セッション作業サマリー

---

**最終更新**: 2024年12月26日  
**ステータス**: Phase 1 月1完了、月2準備完了
