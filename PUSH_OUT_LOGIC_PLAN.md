# 押し出しロジック（Push-Out Logic）実装方針

## 1. 概要

壁面移動により閉じ込められた粒子を、移動する壁面の外側に押し出す処理の実装方針。

## 2. 使用する主要関数

### 2.0 表面衝突モデル（SurfCollide）

```cpp
Particle::OnePart *SurfCollide::collide(
    Particle::OnePart *&ip, 
    double &dtremain, 
    int isurf, 
    double *norm, 
    int isr, 
    int &reaction
)
```

**目的**: 粒子と表面の衝突を処理し、速度を反射モデルに基づいて更新

**パラメータ**:
- `ip`: 粒子へのポインタ（参照渡し、化学反応で削除される可能性あり）
- `dtremain`: 残り時間（押し出し処理ではダミー値として使用）
- `isurf`: 表面要素インデックス
- `norm[3]`: 表面法線ベクトル
- `isr`: 表面反応モデルインデックス（-1の場合は反応なし）
- `reaction`: 出力、実行された反応のインデックス（0の場合は反応なし）

**戻り値**: 化学反応で新しく生成された粒子へのポインタ（通常はNULL）

**実装内容**:
- 各種反射モデル（拡散反射、鏡面反射、CLLモデルなど）を実装
- 表面温度に基づいた速度分布を生成
- 化学反応モデルと連携して反応を処理

**使用箇所**: `move_surf.cpp:Step 5`（押し出し後の速度反射）

**利用可能なモデル**:
- `SurfCollideDiffuse`: 拡散反射（Maxwell分布に基づく）
- `SurfCollideImpulsive`: 衝撃反射（確率的な反射）
- `SurfCollideCLL`: CLLモデル（Cercignani-Lampis-Lordモデル）
- `SurfCollideTD`: 熱脱離モデル
- `SurfCollideAdiabatic`: 断熱反射

### 2.1 表面速度計算

```cpp
void MoveSurf::compute_surf_velocity(int isurf, double *pos, int nevery, double dt, double *vwall)
```

**目的**: 表面要素上の指定位置での速度を計算

**パラメータ**:
- `isurf`: 表面要素インデックス
- `pos[3]`: 表面要素上の位置座標
- `nevery`: 移動ステップ数
- `dt`: タイムステップ
- `vwall[3]`: 出力される速度ベクトル

**実装内容**:
- **並進移動**: `vwall = delta / (nevery * dt)`
- **回転移動**: `vwall = omega × r`（`omega = (theta / dt_total) * rvec`）

**使用箇所**: `move_surf.cpp:1209, 1284`

### 2.2 幾何計算関数

#### 2D: 点から線分への距離

```cpp
double Geometry::distsq_point_line(double *x, double *p1, double *p2)
```

**目的**: 点から線分への最短距離の2乗を計算

**パラメータ**:
- `x[3]`: 粒子位置
- `p1[3]`, `p2[3]`: 線分の端点

**戻り値**: 距離の2乗

**使用箇所**: `move_surf.cpp:1221`

#### 3D: 点から三角形への距離

```cpp
double Geometry::distsq_point_tri(double *x, double *p1, double *p2, double *p3, double *norm)
```

**目的**: 点から三角形への最短距離の2乗を計算

**パラメータ**:
- `x[3]`: 粒子位置
- `p1[3]`, `p2[3]`, `p3[3]`: 三角形の頂点
- `norm[3]`: 三角形の法線ベクトル

**戻り値**: 距離の2乗

**使用箇所**: `move_surf.cpp:1296`

### 2.3 数学ユーティリティ関数

```cpp
// ベクトル減算
void MathExtra::sub3(double *a, double *b, double *c)  // c = a - b

// ベクトル内積
double MathExtra::dot3(double *a, double *b)  // a · b

// ベクトル長さの2乗
double MathExtra::lensq3(double *a)  // |a|²

// 外積
void MathExtra::cross3(double *a, double *b, double *c)  // c = a × b
```

**使用箇所**: 射影パラメータ計算、速度計算など

## 3. 押し出しアルゴリズム（5ステップ）

### Step 1: セル選択（Cell Selection）

**目的**: 粒子が存在するセル内の表面要素を取得

**実装**:
```cpp
int icell = particles[ip].icell;
int nsurf = cells[icell].nsurf;
surfint *csurfs = cells[icell].csurfs;
```

**データ構造**:
- `cells[icell].nsurf`: セル内の表面要素数
- `cells[icell].csurfs`: 表面要素インデックスの配列

**注意点**:
- セルに表面要素がない場合はスキップ（`nno_surf++`）

### Step 2: 速度フィルタ（Velocity Filter）

**目的**: 外側方向に押し出す表面要素のみを選択

**判定条件**:
```cpp
double vn = vwall[0]*norm[0] + vwall[1]*norm[1] + vwall[2]*norm[2];
if (vn <= 0.0) continue;  // 内側方向はスキップ
```

**物理的意味**:
- `vn > 0`: 表面が外側方向に移動 → 粒子を押し出せる
- `vn ≤ 0`: 表面が内側方向に移動 → 粒子を押し出せない

**実装の流れ**:
1. 表面要素の中心点を計算
   - 2D: 線分の中点 `(p1 + p2) / 2`
   - 3D: 三角形の重心 `(p1 + p2 + p3) / 3`
2. `compute_surf_velocity()`で速度を計算
3. 法線ベクトルとの内積で判定

**使用データ**:
- `surf->lines[isurf].norm` (2D)
- `surf->tris[isurf].norm` (3D)

### Step 3: 射影チェック（Projection Check）

**目的**: 粒子から表面要素への最短距離と射影位置を計算

#### 2Dの場合

```cpp
// 距離計算
dist_sq = Geometry::distsq_point_line(x, line->p1, line->p2);

// 射影パラメータ計算（0 ≤ t ≤ 1 が線分上）
double a[3], b[3];
MathExtra::sub3(x, line->p1, a);
MathExtra::sub3(line->p2, line->p1, b);
double len_sq = MathExtra::lensq3(b);
double alpha = MathExtra::dot3(a, b) / len_sq;

// 射影位置
if (alpha >= 0.0 && alpha <= 1.0) {
    // 線分上
    best_pos = p1 + alpha * (p2 - p1);
} else {
    // 端点のうち近い方
    best_pos = (d1_sq < d2_sq) ? p1 : p2;
}
```

#### 3Dの場合

```cpp
// 距離計算
dist_sq = Geometry::distsq_point_tri(x, tri->p1, tri->p2, tri->p3, norm);

// 射影位置（三角形平面への射影）
double a[3];
MathExtra::sub3(x, tri->p1, a);
double alpha = MathExtra::dot3(a, norm);
best_pos = x - alpha * norm;
```

**選択基準**:
- 射影が表面要素の範囲内（`0 ≤ t ≤ 1`）を優先
- 距離が最小の表面要素を選択

### Step 4: 再配置（Re-positioning）

**目的**: 粒子を表面の外側に移動

**実装**:
```cpp
double epsilon = grid->cell_epsilon;  // 安全マージン
if (epsilon <= 0.0) epsilon = 1.0e-9;  // フォールバック値

double new_x[3];
new_x[0] = best_pos[0] + epsilon * best_norm[0];
new_x[1] = best_pos[1] + epsilon * best_norm[1];
new_x[2] = best_pos[2] + epsilon * best_norm[2];

// 位置のみ更新（icellは後で更新）
particles[ip].x[0] = new_x[0];
particles[ip].x[1] = new_x[1];
particles[ip].x[2] = new_x[2];
```

**重要な注意点**:
- `icell`は更新しない（`particle->sort()`で後で更新）
- リンクリスト構造を壊さないため

**安全マージン**:
- `grid->cell_epsilon`: グリッドセルの最小サイズの半分
- 粒子が確実に表面の外側に配置されるように

### Step 5: 速度反射（Velocity Reflection）

**目的**: 粒子の速度を表面衝突モデルに基づいて反射させる

**実装方針**:
1. 表面要素から衝突モデル（`surf_collide`）を取得
2. `surf_collide->collide()`メソッドを呼び出して速度を更新
3. 化学反応モデルが設定されている場合は反応も処理

**実装**:
```cpp
// Step 5: Apply surface collision model to reflect velocity
if (found_candidate && best_surf >= 0) {
    // Get surface element and its collision model
    Surf::Line *line = NULL;
    Surf::Tri *tri = NULL;
    int isc = -1;  // surface collision model index
    double *norm = best_norm;
    
    if (dim == 2) {
        if (best_surf < surf->nlocal) {
            line = &lines[best_surf];
            norm = line->norm;
            isc = line->isc;
        }
    } else {
        if (best_surf < surf->nlocal) {
            tri = &tris[best_surf];
            norm = tri->norm;
            isc = tri->isc;
        }
    }
    
    // Apply collision model if available
    if (isc >= 0 && isc < surf->nsc) {
        SurfCollide *sc = surf->sc[isc];
        if (sc != NULL) {
            // Prepare particle pointer for collision
            Particle::OnePart *ip = &particles[ip];
            
            // Get surface reaction model index (if any)
            int isr = -1;  // surface reaction model index
            if (dim == 2 && line) {
                isr = line->isr;
            } else if (dim == 3 && tri) {
                isr = tri->isr;
            }
            
            // Call collide() method to update velocity
            // Note: dtremain is dummy value since particle is already outside
            double dtremain = update->dt;  // dummy value
            int reaction = 0;
            
            Particle::OnePart *jp = sc->collide(ip, dtremain, best_surf, norm, isr, reaction);
            
            // Handle newly created particle (jp) if any
            // Usually jp is NULL, but if chemistry creates a new particle, handle it
            if (jp != NULL) {
                // New particle created by surface reaction
                // This is rare but should be handled
                // For now, we can ignore or log it
            }
            
            // Update custom properties if needed
            // This is handled by the collide() method internally
        }
    } else {
        // No collision model defined - use specular reflection as fallback
        // v_new = v_old - 2 * (v_old · n) * n
        double *v = particles[ip].v;
        double vn = v[0]*norm[0] + v[1]*norm[1] + v[2]*norm[2];
        v[0] -= 2.0 * vn * norm[0];
        v[1] -= 2.0 * vn * norm[1];
        v[2] -= 2.0 * vn * norm[2];
    }
}
```

**重要な注意点**:
- `surf_collide`モデルが設定されていない場合は、鏡面反射（specular reflection）をフォールバックとして使用
- `collide()`メソッドは粒子の位置が表面に正確にあることを前提としているが、押し出し後は表面の外側にある
- このため、`dtremain`はダミー値として使用（実際の衝突時間計算には使用されない）
- 化学反応モデル（`surf_react`）が設定されている場合、`collide()`内で自動的に処理される

**使用するデータ構造**:
- `surf->lines[isurf].isc`: 2D表面要素の衝突モデルインデックス
- `surf->tris[isurf].isc`: 3D表面要素の衝突モデルインデックス
- `surf->lines[isurf].isr`: 2D表面要素の反応モデルインデックス（-1の場合は反応なし）
- `surf->tris[isurf].isr`: 3D表面要素の反応モデルインデックス（-1の場合は反応なし）
- `surf->sc[isc]`: 衝突モデルへのポインタ

**参考実装箇所**:
- `update.cpp:918-923`: 通常の衝突処理での`surf_collide->collide()`の使用例
- `surf_collide_impulsive.cpp:167-233`: 衝突モデルの実装例
- `surf_collide_diffuse.cpp:140-206`: 拡散反射モデルの実装例
- `surf_collide_cll.cpp:148-239`: CLLモデルの実装例

## 4. 処理フロー

```
1. 検知された粒子リストをループ
   ├─ インデックス検証（安全性チェック）
   ├─ 粒子が存在するセルを取得
   │
2. Step 1: セル内の表面要素を取得
   ├─ nsurf == 0 → スキップ
   │
3. 各表面要素についてループ
   ├─ Step 2: 速度フィルタ
   │  ├─ vn > 0 → 候補
   │  └─ vn ≤ 0 → スキップ
   │
   ├─ Step 3: 射影チェック
   │  ├─ 距離計算
   │  ├─ 射影パラメータ計算
   │  └─ 最適な表面要素を選択
   │
4. Step 4: 再配置
   ├─ 候補が見つかった場合
   │  ├─ 表面位置 + epsilon * 法線
   │  └─ 粒子位置を更新
   └─ 候補が見つからない場合
      └─ nno_candidate++

5. Step 5: 速度反射
   ├─ 候補が見つかった場合
   │  ├─ 表面要素から衝突モデル（isc）を取得
   │  ├─ surf_collide->collide()を呼び出し
   │  │  ├─ 速度を反射モデルに基づいて更新
   │  │  └─ 化学反応があれば処理
   │  └─ 衝突モデルがない場合は鏡面反射を適用
   └─ 候補が見つからない場合
      └─ 速度は変更しない（削除対象）

6. 後処理
   ├─ ゼロ体積セル内の粒子を削除
   └─ 統計情報を集約・出力
```

## 5. エラーハンドリング

### 5.1 安全性チェック

```cpp
// インデックス検証
if (ip < 0 || ip >= particle->nlocal) continue;

// 粒子が削除マークされていないか
if (particles[ip].icell < 0) continue;

// セルインデックスの検証
if (icell < 0 || icell >= grid->nlocal) continue;
```

### 5.2 統計情報

```cpp
int npushed = 0;        // 押し出し成功
int nskipped = 0;       // スキップ（無効インデックスなど）
int nno_surf = 0;       // 表面要素なし
int nno_velocity = 0;   // 速度フィルタで除外
int nno_candidate = 0;  // 候補が見つからない
```

## 6. 後処理

### 6.1 ゼロ体積セルの処理

```cpp
// 押し出し後にゼロ体積セル内の粒子を削除
for (int i = 0; i < nlocal; i++) {
    int icell = particles[i].icell;
    if (icell >= 0 && icell < grid->nlocal) {
        if (cinfo[icell].volume == 0.0) {
            particles[i].icell = -1;  // 削除マーク
            delflag = 1;
        }
    }
}
if (delflag) {
    particle->compress_rebalance();
}
```

### 6.2 セル割り当ての更新

押し出し後、`fix_move_surf.cpp`で以下を実行：

```cpp
// 粒子を新しいセルに割り当て
particle->sort();

// 必要に応じてプロセス間移動
comm->migrate_particles();
```

## 7. 実装上の注意点

### 7.1 メモリ安全性

- `trapped_list`のインデックスは`compress_rebalance()`後に無効になる可能性がある
- ループ内で`particles`ポインタを再取得する

### 7.2 並列計算

- 各プロセスが自分の所有粒子のみ処理
- MPI集約で統計情報を集計

### 7.3 数値精度

- `epsilon`は十分小さく、かつ粒子が確実に外側に配置される値
- 距離計算は2乗で比較（平方根を避ける）

## 8. 関数設計の詳細

### 8.1 `push_out_particles()`関数の拡張

**現在の実装**: `move_surf.cpp:1113-1384`

**拡張内容**:
- Step 4（再配置）の後にStep 5（速度反射）を追加
- 表面要素から衝突モデルを取得
- `surf_collide->collide()`を呼び出して速度を更新

**関数シグネチャ**（変更なし）:
```cpp
void MoveSurf::push_out_particles(int *trapped_list, int ntrapped, int nevery, double dt)
```

**追加する処理**:
```cpp
// Step 5: Velocity reflection (after Step 4)
if (found_candidate && best_surf >= 0) {
    apply_velocity_reflection(&particles[ip], best_surf, best_norm);
}
```

### 8.2 補助関数: `apply_velocity_reflection()`

**目的**: 粒子の速度を表面衝突モデルに基づいて反射させる

**シグネチャ**:
```cpp
void MoveSurf::apply_velocity_reflection(
    Particle::OnePart *p, 
    int isurf, 
    double *norm
)
```

**実装方針**:
1. 表面要素の種類（2D/3D）を判定
2. 衝突モデルインデックス（`isc`）と反応モデルインデックス（`isr`）を取得
3. `surf_collide->collide()`を呼び出し
4. 衝突モデルがない場合は鏡面反射を適用

**実装場所**: `move_surf.cpp`内に新規追加

**注意点**:
- `collide()`メソッドは粒子が表面に正確にあることを前提としている
- 押し出し後は粒子が表面の外側にあるため、`dtremain`はダミー値として使用
- 化学反応で粒子が削除される可能性があるが、押し出し処理では通常発生しない

### 8.3 実装の統合ポイント

**`move_surf.cpp`内の変更箇所**:
- `push_out_particles()`関数内、Step 4の後にStep 5を追加
- 新しい補助関数`apply_velocity_reflection()`を追加

**`move_surf.h`内の変更箇所**:
- `apply_velocity_reflection()`の宣言を追加（プライベートメソッド）

**必要なインクルード**:
```cpp
#include "surf_collide.h"
#include "surf_react.h"
```

## 9. 参考実装箇所

- **押し出し関数定義**: `move_surf.cpp:1113-1384`
- **呼び出し**: `fix_move_surf.cpp:270`
- **速度計算**: `move_surf.cpp:1063-1101`
- **幾何計算**: `geometry.h`内の`Geometry`クラス
- **衝突処理の参考**: `update.cpp:918-923`（通常の衝突処理）
- **衝突モデル実装**: `surf_collide_impulsive.cpp:167-233`
- **拡散反射モデル**: `surf_collide_diffuse.cpp:140-206`
- **CLLモデル**: `surf_collide_cll.cpp:148-239`

## 10. 実装の詳細ロジック

### 10.1 Step 5の実装フロー

```
1. 最適な表面要素（best_surf）が決定された後
   │
2. 表面要素の種類を判定（2D/3D）
   ├─ 2D: line = &lines[best_surf]
   │  ├─ isc = line->isc
   │  └─ isr = line->isr
   └─ 3D: tri = &tris[best_surf]
      ├─ isc = tri->isc
      └─ isr = tri->isr
   │
3. 衝突モデルの存在確認
   ├─ isc >= 0 && isc < surf->nsc → 衝突モデルあり
   └─ それ以外 → 衝突モデルなし（鏡面反射を使用）
   │
4. 衝突モデルがある場合
   ├─ sc = surf->sc[isc]
   ├─ dtremain = update->dt（ダミー値）
   ├─ sc->collide(ip, dtremain, best_surf, norm, isr, reaction)
   │  ├─ 速度を反射モデルに基づいて更新
   │  ├─ 化学反応があれば処理（isr >= 0の場合）
   │  └─ 新粒子（jp）が生成された場合は処理
   └─ 完了
   │
5. 衝突モデルがない場合（フォールバック）
   ├─ 鏡面反射を適用
   ├─ v_new = v_old - 2 * (v_old · n) * n
   └─ 完了
```

### 10.2 エッジケースの処理

**ケース1: 衝突モデルが設定されていない**
- 解決策: 鏡面反射（specular reflection）をフォールバックとして使用
- 実装: 法線方向の速度成分を反転

**ケース2: 化学反応で粒子が削除される**
- 解決策: `collide()`メソッドが`ip`をNULLに設定する可能性がある
- 実装: 押し出し処理では通常発生しないが、チェックを追加
```cpp
if (ip == NULL) {
    // Particle was destroyed by chemistry
    // Mark for deletion (already handled by collide())
    continue;
}
```

**ケース3: 化学反応で新粒子が生成される**
- 解決策: `collide()`メソッドが`jp`を返す可能性がある
- 実装: 通常はNULLだが、チェックを追加してログ出力
```cpp
if (jp != NULL) {
    // New particle created - log for debugging
    // Usually this doesn't happen in push-out context
}
```

**ケース4: 表面要素のインデックスが無効**
- 解決策: インデックス範囲チェックを追加
```cpp
if (best_surf < 0 || best_surf >= surf->nlocal) {
    // Invalid surface index - skip reflection
    continue;
}
```

### 10.3 デバッグとログ出力

**追加する統計情報**:
```cpp
int nreflected = 0;        // 速度反射成功
int nspecular = 0;         // 鏡面反射（フォールバック）
int nno_collide = 0;       // 衝突モデルなし
int nchemistry = 0;        // 化学反応が発生
```

**ログ出力例**:
```cpp
if (me == 0) {
    fprintf(screen, "Reflection stats: reflected=%ld, specular=%ld, no_collide=%ld, chemistry=%ld\n",
            nreflected_global, nspecular_global, nno_collide_global, nchemistry_global);
}
```

## 11. 今後の改善点（予定）

1. **速度計算の精度向上**
   - 表面要素上の実際の位置での速度計算（現在は中心点）
   - 回転の場合、位置依存の速度を正確に計算

2. **射影アルゴリズムの改善**
   - 3Dでの三角形内判定の精度向上
   - エッジケースの処理改善

3. **複数候補の優先順位**
   - 距離だけでなく、速度成分も考慮
   - より物理的に適切な表面要素を選択

4. **エラーハンドリングの強化**
   - 候補が見つからない場合の代替処理
   - デバッグ情報の充実

5. **速度反射の改善**
   - 移動する壁面との相対速度を考慮した反射
   - より物理的に正確な反射モデルの適用

---

**最終更新**: 2024年12月26日  
**ステータス**: 実装方針確定（速度反射処理追加）

## 12. 実装チェックリスト

### 12.1 必須実装項目

- [ ] Step 5（速度反射）の実装
  - [ ] `apply_velocity_reflection()`補助関数の実装
  - [ ] `push_out_particles()`内でのStep 5呼び出し
  - [ ] 衝突モデルの取得と使用
  - [ ] 鏡面反射フォールバックの実装
  - [ ] エッジケースの処理

- [ ] ヘッダーファイルの更新
  - [ ] `move_surf.h`に`apply_velocity_reflection()`の宣言を追加
  - [ ] 必要なインクルードの追加

- [ ] 統計情報の追加
  - [ ] 反射統計のカウント
  - [ ] ログ出力の追加

### 12.2 テスト項目

- [ ] 衝突モデルが設定されている場合の動作確認
- [ ] 衝突モデルが設定されていない場合の鏡面反射動作確認
- [ ] 化学反応モデルが設定されている場合の動作確認
- [ ] 複数の表面要素がある場合の最適選択の確認
- [ ] 並列計算環境での動作確認
