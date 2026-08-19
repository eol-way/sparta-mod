# SPARTA改造計画書 - 実装詳細版

**作成日**: 2026年01月30日  
**参照元**: `build/SPARTA改造計画書Rev1.pdf`  
**目的**: 壁内粒子の押し出し処理と回転壁の速度考慮を実装

---

## 📋 目次

1. [プロジェクト概要](#プロジェクト概要)
2. [開発方針とルール](#開発方針とルール)
3. [Phase 1: 粒子押し出し (Sweeping)](#phase-1-粒子押し出し-sweeping)
4. [Phase 2: 壁速度補正の運用整理](#phase-2-壁速度補正の運用整理)
5. [実装チェックリスト](#実装チェックリスト)
6. [デバッグとテスト](#デバッグとテスト)

---

## プロジェクト概要

### 目的
壁移動時に壁内に取り残された粒子を削除するのではなく、**幾何学的に最近傍の壁面へ射影して押し出し、反射速度を更新する**機能を実装する。

### 現状の問題点
- 現在の`MoveSurf::remove_particles()`は壁内粒子を単純に削除している
- 回転壁の速度が反射時に考慮されていない

### 変更の影響範囲
- **主要ファイル**: `src/move_surf.cpp`, `src/move_surf.h`
- **関連ファイル**: `src/surf_collide_*.cpp`, `src/update.cpp`, `src/grid_surf.cpp`
- **新規ファイル**: なし（既存ファイルの拡張）

---

## 開発方針とルール

### 1. コード変更のトレーサビリティ

**統一タグルール**:
- すべての変更箇所に `// [MODIFIED] Phase X` タグを付与
- 削除したコードはコメントアウトで残し、`// [REMOVED] Phase X: 理由` を記載
- 新規追加コードには `// [ADDED] Phase X: 説明` を記載

**例**:
```cpp
// [MODIFIED] Phase 1: セル内無差別削除判定から幾何学判定に変更
if (geometric_distance_check(...)) {
  // [ADDED] Phase 1: 最近傍点への射影処理
  project_to_nearest_surface(...);
}
```

### 2. デバッグログの埋め込み

**デバッグフラグ**:
```cpp
#define DEBUG_SPARTA_SWEEP 0  // Phase 1用
#define DEBUG_SPARTA_ROTATE 0  // Phase 2用
```

**ログ出力例**:
```cpp
#if DEBUG_SPARTA_SWEEP
  if (me == 0) {
    fprintf(screen, "SWEEP: Particle %d at (%g,%g,%g) projected to (%g,%g,%g)\n",
            particle_id, x_old[0], x_old[1], x_old[2], x_new[0], x_new[1], x_new[2]);
  }
#endif
```

### 3. 既存コードの保持

- **削除せずコメントアウト**: 既存の`remove_particles()`ロジックは残す
- **後方互換性**: 既存の動作モードを維持（フラグで切り替え可能に）

---

## Phase 1: 粒子押し出し (Sweeping)

### 1.1 概要

壁移動後に壁内に取り残された粒子を検知し、最近傍の非透明な壁面要素へ射影して押し出す。純粋な幾何学的判定のみを使用。透明な表面要素は除外（粒子は通過するため、壁内判定・反射処理の対象外）。

### 1.2 実装タスク

#### [ ] タスク 1.1: 壁内粒子の検知機能の実装

**ファイル**: `src/move_surf.cpp`, `src/move_surf.h`

**関数名**: `detect_inside_particles()` (新規)

**実装内容**:
1. INSIDE/OVERLAPセルを先に特定し、それらのセル内の粒子のみを走査する（OUTSIDEセルの粒子は触らないため高速化）
2. セルごとに`cinfo[icell].first`と`particle->next[i]`のリンクリストで粒子をたどる
3. INSIDEセルの粒子は確実に壁内として検出
4. OVERLAPセル（非透明な表面要素を持つセル）の粒子は`Grid::outside_surfs()`を使用して壁内判定
5. 透明な表面要素のみのセルは除外（粒子は通過するため）
6. 子セル（`nsplit <= 0`）は除外し、親/unsplitセルのみを探索
7. 壁内粒子のリストを作成

**注意**: `remove_particles()`呼び出し前に`particle->sort()`が実行されていること（`command()`内で保証）を前提とする。

**関数シグネチャ**:
```cpp
// src/move_surf.h
private:
  struct InsideParticle {
    int particle_index;
    int icell;
    double x[3];
  };
  std::vector<InsideParticle> detect_inside_particles();
```

**実装詳細**:
```cpp
// src/move_surf.cpp
std::vector<MoveSurf::InsideParticle> MoveSurf::detect_inside_particles()
{
  // [ADDED] Phase 1: 壁内粒子の検知（セル先走査でOUTSIDEセル内粒子をスキップ）
  std::vector<InsideParticle> inside_particles;

  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Particle::OnePart *particles = particle->particles;
  int *next = particle->next;
  int nglocal = grid->nlocal;

  for (int icell = 0; icell < nglocal; icell++) {
    if (cells[icell].nsplit <= 0) continue;  // 子セルは除外（親/unsplitのみ探索）
    if (cinfo[icell].type != INSIDE && cinfo[icell].type != OVERLAP) continue;

    int i = cinfo[icell].first;
    while (i >= 0) {
      if (cinfo[icell].type == INSIDE) {
        InsideParticle p;
        p.particle_index = i;
        p.icell = icell;
        p.x[0] = particles[i].x[0];
        p.x[1] = particles[i].x[1];
        p.x[2] = particles[i].x[2];
        inside_particles.push_back(p);
      } else {
        // OVERLAPセル: 幾何学的判定
        double xcell[3];
        if (grid->point_outside_surfs(icell, xcell)) {
          if (!grid->outside_surfs(icell, particles[i].x, xcell)) {
            InsideParticle p;
            p.particle_index = i;
            p.icell = icell;
            p.x[0] = particles[i].x[0];
            p.x[1] = particles[i].x[1];
            p.x[2] = particles[i].x[2];
            inside_particles.push_back(p);
          }
        }
      }
      i = next[i];
    }
  }

  return inside_particles;
}
```

---

#### [ ] タスク 1.2: 最近傍点への射影処理の実装

**ファイル**: `src/move_surf.cpp`, `src/geometry.cpp` (既存関数を利用)

**関数名**: `project_to_nearest_surface()` (新規)

**実装内容**:
1. 粒子位置から最近傍の非透明な表面要素（Line/Tri）を探索
2. 透明な表面要素はスキップ（粒子は通過するため）
3. 表面要素上の最近傍点を計算
4. 安全マージンεを付与して再配置

**関数シグネチャ**:
```cpp
// src/move_surf.h
private:
  int project_to_nearest_surface(double *x, int icell, double *x_projected,
                                 double *norm, int &isurf);
  double epsilon_sweep;  // 安全マージン（デフォルト: 1e-6）
```

**実装詳細**:
```cpp
// src/move_surf.cpp
int MoveSurf::project_to_nearest_surface(double *x, int icell,
                                         double *x_projected,
                                         double *norm, int &isurf)
{
  // [ADDED] Phase 1: 最近傍表面への射影
  
  double mindist = BIG;
  int minisurf = -1;
  double minproj[3], minnorm[3];

  
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;

  // 呼び出し側から受け取ったセルIDを使用（xから逆引きしない）
  if (icell < 0 || icell >= grid->nlocal) return 0;

  // 子セルなら親splitセルへ正規化してから探索（なぜなら、nsurf = cells[icell].nsurfで呼び出されるのは親セルが所有するサーフェスになるから）
  if (cells[icell].nsplit <= 0)
    icell = sinfo[cells[icell].isplit].icell;
  
  // OVERLAPセルのみ処理（非透明な表面要素を持つセル）
  if (cinfo[icell].type == OVERLAP && cells[icell].nsurf > 0) {
    surfint *csurfs = cells[icell].csurfs;
    int nsurf = cells[icell].nsurf;
    
    for (int m = 0; m < nsurf; m++) {
      int isurf_local = csurfs[m];
      
      // 透明な表面要素はスキップ（粒子は通過するため）
      if (dim == 2) {
        if (surf->lines[isurf_local].transparent) continue;
      } else {
        if (surf->tris[isurf_local].transparent) continue;
      }

      double dist;
      if (dim == 2) {
        Surf::Line *line = &surf->lines[isurf_local];
        dist = Geometry::distsq_point_line(x, line->p1, line->p2);
      } else {
        Surf::Tri *tri = &surf->tris[isurf_local];
        dist = Geometry::distsq_point_tri(x, tri->p1, tri->p2, tri->p3, tri->norm);
      }

      // ループ内は最小距離探索のみ
      if (dist < mindist) {
        mindist = dist;
        minisurf = isurf_local;
      }
    }
  }
  
  if (minisurf >= 0) {
    if (dim == 2) {
      Surf::Line *line = &surf->lines[minisurf];

      // 法線は事前計算済み値を再利用（毎回の外積計算はしない）
      minnorm[0] = line->norm[0];
      minnorm[1] = line->norm[1];
      minnorm[2] = line->norm[2];

      // 点xを線分(p1-p2)へ射影し、tを[0,1]にクランプして線分上にスナップ
      double dx = line->p2[0] - line->p1[0];
      double dy = line->p2[1] - line->p1[1];
      double dz = line->p2[2] - line->p1[2];
      double len2 = dx*dx + dy*dy + dz*dz;

      if (len2 > 0.0) {
        double t = ((x[0] - line->p1[0]) * dx +
                    (x[1] - line->p1[1]) * dy +
                    (x[2] - line->p1[2]) * dz) / len2;
        if (t < 0.0) t = 0.0;
        else if (t > 1.0) t = 1.0;

        minproj[0] = line->p1[0] + t * dx;
        minproj[1] = line->p1[1] + t * dy;
        minproj[2] = line->p1[2] + t * dz;
      } else {
        // 退化線分（長さ0）は端点へ退避
        minproj[0] = line->p1[0];
        minproj[1] = line->p1[1];
        minproj[2] = line->p1[2];
      }
    } else {
      Surf::Tri *tri = &surf->tris[minisurf];

      // 法線は事前計算済み値を再利用
      minnorm[0] = tri->norm[0];
      minnorm[1] = tri->norm[1];
      minnorm[2] = tri->norm[2];

      // 点xを三角形が張る無限平面へ射影（垂線の足）
      double vx = x[0] - tri->p1[0];
      double vy = x[1] - tri->p1[1];
      double vz = x[2] - tri->p1[2];
      double signed_dist = vx*minnorm[0] + vy*minnorm[1] + vz*minnorm[2];

      minproj[0] = x[0] - signed_dist * minnorm[0];
      minproj[1] = x[1] - signed_dist * minnorm[1];
      minproj[2] = x[2] - signed_dist * minnorm[2];
    }

    x_projected[0] = minproj[0];
    x_projected[1] = minproj[1];
    x_projected[2] = minproj[2];
    norm[0] = minnorm[0];
    norm[1] = minnorm[1];
    norm[2] = minnorm[2];

    // 安全マージンεを法線方向に付与
    x_projected[0] += epsilon_sweep * norm[0];
    x_projected[1] += epsilon_sweep * norm[1];
    x_projected[2] += epsilon_sweep * norm[2];
    isurf = minisurf;
    return 1;
  }
  
  return 0;  // 射影失敗
}
```

**使用する既存関数**:
- `distsq_point_line()` (geometry.cpp:1251)
- `distsq_point_tri()` (geometry.cpp:1289)

---

#### [ ] タスク 1.3: 反射速度の更新処理

**ファイル**: `src/move_surf.cpp`, `src/surf_collide_*.cpp`

**関数名**: `update_reflection_velocity()` (新規)

**実装内容**:
1. 射影した表面要素の`isc`（surface collide ID）を取得
2. `surf->sc[isc]` の仮想関数 `collide()` を呼び出し
3. 反射速度を計算して粒子速度を更新（モデル選択はSPARTAの仮想関数に委譲）

**関数シグネチャ**:
```cpp
// src/move_surf.h
private:
  void update_reflection_velocity(Particle::OnePart *p, int isurf, double *norm);
```

**実装詳細**:
```cpp
// src/move_surf.cpp
void MoveSurf::update_reflection_velocity(Particle::OnePart *p, int isurf,
                                          double *norm)
{
  // [ADDED] Phase 1: 反射速度の更新
  
  Surf::Line *lines = surf->lines;
  Surf::Tri *tris = surf->tris;
  
  int isc;  // surface collide ID
  if (dim == 2) {
    isc = lines[isurf].isc;
  } else {
    isc = tris[isurf].isc;
  }

  // セーフティチェック: 反射モデル未割当て(-1)や範囲外アクセスを防ぐ
  if (isc < 0 || isc >= surf->nsc || surf->sc[isc] == NULL) {
    if (me == 0) {
      error->warning(FLERR,
        "Surface collide model is not assigned for projected surface element");
    }
    return;
  }

  // SPARTAの仮想関数でモデル選択を委譲（Diffuse/Specular/CLLをif分岐しない）
  // surf_modifyで割り当てられたiscに応じて適切な派生クラスのcollide()が呼ばれる
  double dt_dummy = 0.0;
  int isr = -1;  // 化学反応なし
  int reaction = 0;
  surf->sc[isc]->collide(p, dt_dummy, isurf, norm, isr, reaction);
  //現状のpには反射前粒子速度がそのまま入っているので、specular処理でも問題ない
  
  // Phase 2方針:
  // 壁速度補正は surf->sc[isc]->collide() の内部実装
  // （surf_collide diffuse/cll の translate/rotate 等）に委譲する。
  // MoveSurf 側で追加の速度加算は行わない。
}
```

**注意事項**:
- `SurfCollide` 派生クラス（Diffuse/Specular/CLLなど）の直接メソッド呼び出しは行わない
- `surf->sc[isc]->collide()` の仮想関数ディスパッチにモデル選択を委譲する
- `isc` は未割当てで `-1` の可能性があるため、呼び出し前に範囲・NULLチェックを行う
- `isr = -1`（化学反応なし）では戻り値 `jp` は使用しないため、捨て打ちで呼び出す
- `collide()` 呼び出しにより `nsingle` など統計カウンタは増加する（コストは軽微）


---

#### [ ] タスク 1.4: INSIDEセル粒子の救済処理

**ファイル**: `src/move_surf.cpp`

**関数名**: `rescue_inside_particles()` (新規)

**実装内容**:
1. INSIDEセルの粒子に対して近傍探索を実行
2. 近傍セルリング（Ring 0→1→2）内のOVERLAPセルを探索
3. 透明な表面要素はスキップ（粒子は通過するため、救済対象外）
4. 救済可能な場合は射影、不可能な場合は削除マークして警告ログ

**関数シグネチャ**:
```cpp
// src/move_surf.h
private:
  int rescue_inside_particles(std::vector<InsideParticle> &inside_particles,
                              int &ndeleted_marked);
```

**実装詳細**:
```cpp
// src/move_surf.cpp
int MoveSurf::rescue_inside_particles(std::vector<InsideParticle> &inside_particles,
                                      int &ndeleted_marked)
{
  // [ADDED] Phase 1: INSIDEセル粒子の救済
  
  int nrescued = 0;
  ndeleted_marked = 0;
  
  for (auto &p : inside_particles) {
    int icell = p.icell;

    double x_proj[3], norm[3];
    int isurf;

    // 局所リング探索（簡易・高速）:
    //   Ring 0 = 現在セル, Ring 1 = 面隣接, Ring 2 = さらに外側
    //   OVERLAPかつ非透明surfaceのみ候補にする
    const int max_ring = 2;
    if (find_nearest_surface_in_rings(p.x, icell, max_ring, x_proj, norm, isurf)) {
      // 救済成功: 射影して速度更新
      Particle::OnePart *part = &particle->particles[p.particle_index];
      part->x[0] = x_proj[0];
      part->x[1] = x_proj[1];
      part->x[2] = x_proj[2];
      
      update_reflection_velocity(part, isurf, norm);
      nrescued++;
      
#if DEBUG_SPARTA_SWEEP
      if (me == 0) {
        fprintf(screen, "RESCUE: Particle %d rescued from INSIDE cell %d\n",
                p.particle_index, icell);
      }
#endif
    } else {
      // 救済失敗: 削除マークして警告
      Particle::OnePart *part = &particle->particles[p.particle_index];
      part->icell = -1;  // compress_rebalance()で一括削除
      ndeleted_marked++;
      
      if (me == 0) {
        error->warning(FLERR, 
          "Particle in INSIDE cell could not be rescued, deleted");
      }
    }
  }
  
  return nrescued;
}
```

**補助関数（追加）**:
```cpp
// src/move_surf.cpp
int MoveSurf::find_nearest_surface_in_rings(double *x, int icell, int max_ring,
                                            double *x_proj, double *norm, int &isurf)
{
  Grid::ChildCell *cells = grid->cells;
  Grid::ChildInfo *cinfo = grid->cinfo;
  Grid::SplitInfo *sinfo = grid->sinfo;
  Grid::ParentCell *pcells = grid->pcells;

  const int nall = grid->nlocal + grid->nghost;
  if (icell < 0 || icell >= nall) return 0;

  // subcell -> 親splitセルへ正規化
  if (cells[icell].nsplit <= 0) {
    int isplit = cells[icell].isplit;
    if (isplit < 0) return 0;
    icell = sinfo[isplit].icell;
    if (icell < 0 || icell >= nall) return 0;
  }

  // 訪問管理は「訪問したセルIDのみ」を小さく保持
  std::vector<int> visited_list;
  visited_list.reserve(64);

  auto is_visited = [&visited_list](int c) {
    for (int v : visited_list) if (v == c) return 1;
    return 0;
  };

  std::vector<int> frontier, next_frontier;
  frontier.reserve(32);
  next_frontier.reserve(64);
  frontier.push_back(icell);
  visited_list.push_back(icell);

  const int nface = (dim == 3) ? 6 : 4;
  const double BIGDIST = 1.0e300;
  const double eps = 0.25 * grid->cell_epsilon;

  for (int ring = 0; ring <= max_ring && !frontier.empty(); ring++) {
    double best_d2 = BIGDIST;
    int best_isurf = -1;
    double best_proj[3] = {0.0, 0.0, 0.0};
    double best_norm[3] = {0.0, 0.0, 0.0};

    // このRing内で候補表面を評価
    for (int c : frontier) {
      // cinfoはownedセルにしかないためguardする
      if (c < 0 || c >= grid->nlocal) continue;
      if (cinfo[c].type != OVERLAP) continue;
      if (cells[c].nsurf <= 0) continue;

      double cand_proj[3], cand_norm[3];
      int cand_isurf = -1;
      if (project_to_nearest_surface(x, c, cand_proj, cand_norm, cand_isurf)) {
        double dx = x[0] - cand_proj[0];
        double dy = x[1] - cand_proj[1];
        double dz = x[2] - cand_proj[2];
        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < best_d2) {
          best_d2 = d2;
          best_isurf = cand_isurf;
          best_proj[0] = cand_proj[0];
          best_proj[1] = cand_proj[1];
          best_proj[2] = cand_proj[2];
          best_norm[0] = cand_norm[0];
          best_norm[1] = cand_norm[1];
          best_norm[2] = cand_norm[2];
        }
      }
    }

    // 近場優先: このRingで見つかったら即採用
    if (best_isurf >= 0) {
      isurf = best_isurf;
      x_proj[0] = best_proj[0];
      x_proj[1] = best_proj[1];
      x_proj[2] = best_proj[2];
      norm[0] = best_norm[0];
      norm[1] = best_norm[1];
      norm[2] = best_norm[2];
      return 1;
    }

    // 次Ringを構築（面隣接）
    next_frontier.clear();
    for (int c : frontier) {
      for (int iface = 0; iface < nface; iface++) {
        int nflag = grid->neigh_decode(cells[c].nmask, iface);
        int j = -1;

        if (nflag == 0 || nflag == 3) {
          // NCHILD / NPBCHILD: neigh[] は cells index
          j = static_cast<int>(cells[c].neigh[iface]);

        } else if (nflag == 1 || nflag == 4) {
          // NPARENT / NPBPARENT: neigh[] は pcells index
          int ip = static_cast<int>(cells[c].neigh[iface]);
          if (ip < 0 || ip >= grid->nparent) continue;
          Grid::ParentCell *pc = &pcells[ip];

          // face外側へ少し出したprobe点で子セルindexを解決
          double probe[3];
          probe[0] = 0.5 * (cells[c].lo[0] + cells[c].hi[0]);
          probe[1] = 0.5 * (cells[c].lo[1] + cells[c].hi[1]);
          probe[2] = 0.5 * (cells[c].lo[2] + cells[c].hi[2]);
          if (iface == 0) probe[0] = cells[c].lo[0] - eps;
          else if (iface == 1) probe[0] = cells[c].hi[0] + eps;
          else if (iface == 2) probe[1] = cells[c].lo[1] - eps;
          else if (iface == 3) probe[1] = cells[c].hi[1] + eps;
          else if (iface == 4) probe[2] = cells[c].lo[2] - eps;
          else if (iface == 5) probe[2] = cells[c].hi[2] + eps;

          j = grid->id_find_child(pc->id, cells[c].level, pc->lo, pc->hi, probe);
        } else {
          // NUNKNOWN / NPBUNKNOWN / NBOUND は探索しない
          continue;
        }

        if (j < 0 || j >= nall) continue;

        // subcellなら親splitセルへ正規化
        if (cells[j].nsplit <= 0) {
          int jsplit = cells[j].isplit;
          if (jsplit < 0) continue;
          j = sinfo[jsplit].icell;
          if (j < 0 || j >= nall) continue;
        }

        if (!is_visited(j)) {
          visited_list.push_back(j);
          next_frontier.push_back(j);
        }
      }
    }

    frontier.swap(next_frontier);
  }

  return 0;
}
```

---

#### [ ] タスク 1.5: `remove_particles()`の改修

**ファイル**: `src/move_surf.cpp`

**変更内容**:
1. 既存の削除ロジックをコメントアウト（後方互換性のため）
2. 新しい押し出し処理を実装
3. フラグで切り替え可能に

**実装詳細**:
```cpp
// src/move_surf.cpp:819-885
bigint MoveSurf::remove_particles()
{
  // [MODIFIED] Phase 1: 削除から押し出し処理に変更
  
  // [REMOVED] Phase 1: 既存の削除ロジック（コメントアウトで保持）
  /*
  int isurf,nsurf;
  surfint *csurfs;
  // ... 既存コード ...
  */
  
  // [ADDED] Phase 1: 新しい押し出し処理
  int delflag = 0;
  int nlocal_old = particle->nlocal;
  int sweep_mode = 1;  // 0=削除モード（旧）, 1=押し出しモード（新）
  
  if (sweep_mode == 0) {
    // 旧モード: 削除処理（既存コードを実行）
    // ...
  } else {
    // 新モード: 押し出し処理
    std::vector<InsideParticle> inside_particles = detect_inside_particles();
    std::vector<InsideParticle> rescue_targets;
    rescue_targets.reserve(inside_particles.size());
    
    int nswept = 0;
    int nrescued = 0;
    int ndeleted = 0;
    
    for (auto &p : inside_particles) {
      Particle::OnePart *part = &particle->particles[p.particle_index];
      
      double x_proj[3], norm[3];
      int isurf;
      
      if (project_to_nearest_surface(part->x, p.icell, x_proj, norm, isurf)) {
        // 射影成功: 位置と速度を更新
        part->x[0] = x_proj[0];
        part->x[1] = x_proj[1];
        part->x[2] = x_proj[2];
        
        update_reflection_velocity(part, isurf, norm);
        nswept++;
      } else {
        // 射影失敗: 一旦キューへ積んでまとめて救済
        rescue_targets.push_back(p);
      }
    }

    if (!rescue_targets.empty()) {
      int ndeleted_marked = 0;
      nrescued = rescue_inside_particles(rescue_targets, ndeleted_marked);
      ndeleted += ndeleted_marked;
      if (ndeleted_marked > 0) delflag = 1;
    }
    
    // 統計情報を出力
    if (me == 0) {
      if (screen) {
        fprintf(screen, "SWEEP: %d particles swept, %d rescued, %d deleted\n",
                nswept, nrescued, ndeleted);
      }
    }
  }
  
  // 粒子の再圧縮
  if (delflag) particle->compress_rebalance();
  
  // 削除数を返す（後方互換性）
  bigint delta = nlocal_old - particle->nlocal;
  bigint ndeleted_total;
  MPI_Allreduce(&delta, &ndeleted_total, 1, MPI_SPARTA_BIGINT, MPI_SUM, world);
  return ndeleted_total;
}
```

---

#### [ ] タスク 1.6: ヘッダーファイルの更新

**ファイル**: `src/move_surf.h`

**追加内容**:
1. 新しい構造体・関数の宣言
2. メンバ変数の追加

**実装詳細**:
```cpp
// src/move_surf.h
#include <vector>
#include "particle.h"

class MoveSurf : protected Pointers {
  // ... 既存のpublicメンバ ...
  
 private:
  // ... 既存のprivateメンバ ...
  
  // [ADDED] Phase 1: 押し出し処理用のメンバ
  struct InsideParticle {
    int particle_index;
    int icell;
    double x[3];
  };
  
  double epsilon_sweep;  // 安全マージン（デフォルト: 1e-6）
  
  std::vector<InsideParticle> detect_inside_particles();
  int project_to_nearest_surface(double *x, int icell, double *x_projected,
                                 double *norm, int &isurf);
  void update_reflection_velocity(Particle::OnePart *p, int isurf, double *norm);
  int rescue_inside_particles(std::vector<InsideParticle> &inside_particles,
                              int &ndeleted_marked);
  int find_nearest_surface_in_rings(double *x, int icell, int max_ring,
                                    double *x_proj, double *norm, int &isurf);
};
```

**注意事項**:
- `find_nearest_surface_in_rings()` は透明な表面要素をスキップする必要がある
- 探索対象はOVERLAPセル内の非透明な表面要素のみとする
- `icell` が子セルの場合は親splitセルへ正規化してからリング探索を開始する
- 探索範囲は固定上限（例: `max_ring = 2`）で打ち切り、計算コストを抑える
- `std::vector` と `Particle::OnePart` を使うため、`move_surf.h` に `#include <vector>` と `#include "particle.h"` を追加する
- `epsilon_sweep` はコンストラクタで初期値（例: `1.0e-6`）を設定する

---

## Phase 2: 壁速度補正の運用整理

### 2.1 概要

押し出し反射時の速度補正は、`update_reflection_velocity()` 内で呼ぶ
`surf->sc[isc]->collide()` の既存ロジックに委譲する。  
`MoveSurf`/`FixMoveSurf` 側で独自の壁速度（角速度・並進速度）を算出して加算する処理は導入しない。

### 2.2 実装タスク

#### [ ] タスク 2.1: `surf_collide` 側補正の利用方針を固定

**ファイル**: `src/move_surf.cpp`, `src/surf_collide_*.cpp`

**実装内容**:
1. `update_reflection_velocity()` では `surf->sc[isc]->collide()` のみを呼び、追加加算しない
2. 壁速度補正は `surf_collide` モデル側（例: diffuse/cll の `translate/rotate`）に委譲する
3. 二重加算回避のため、MoveSurf/FixMoveSurf 側に独自補正コードを追加しない

---

#### [ ] タスク 2.2: `update_reflection_velocity()` の責務を固定

**ファイル**: `src/move_surf.cpp`

**変更箇所**: `update_reflection_velocity()`

**実装内容**:
1. `surf->sc[isc]->collide()` を唯一の反射・補正入口として使用
2. MoveSurf 側での `p->v` 追加加算を行わない
3. 反射モデル切替は仮想関数ディスパッチに委譲

**実装詳細**:
```cpp
// src/move_surf.cpp
void MoveSurf::update_reflection_velocity(Particle::OnePart *p, int isurf,
                                          double *norm)
{
  // ... Phase 1の実装と同一 ...
  // surf->sc[isc]->collide() が速度補正を含む反射を実施する。
}
```

---

#### [ ] タスク 2.3: `fix move/surf` での運用確認

**ファイル**: `src/fix_move_surf.cpp`

**変更内容**:
1. `fix move/surf` の押し出し反射でも `update_reflection_velocity()` を経由する
2. 追加の kinematics 受け渡しは実装しない
3. 必要な補正は `surf_collide` 入力設定（`translate/rotate`）側で管理する

**実装詳細**:
```cpp
// src/fix_move_surf.cpp
// [MODIFIED] Phase 2: fix move/surfでも同じ反射経路を使用
// end_of_step() -> movesurf->remove_particles()
// -> update_reflection_velocity() -> surf->sc[isc]->collide()
```

**注意**:
- 押し出し反射で壁速度補正を使う場合は、対象 `surf_collide` モデルに
  `translate/rotate` を設定する。
- MoveSurf/FixMoveSurf 側には独自速度補正を実装しない（責務分離）。

---

## 実装チェックリスト

### Phase 1: 粒子押し出し

- [ ] **タスク 1.1**: 壁内粒子の検知機能の実装
- [ ] **タスク 1.2**: 最近傍点への射影処理の実装
- [ ] **タスク 1.3**: 反射速度の更新処理
- [ ] **タスク 1.4**: INSIDEセル粒子の救済処理
- [ ] **タスク 1.5**: `remove_particles()`の改修
- [ ] **タスク 1.6**: ヘッダーファイルの更新

### Phase 2: 壁速度補正の運用整理

- [ ] **タスク 2.1**: `surf_collide` 側補正の利用方針を固定
- [ ] **タスク 2.2**: `update_reflection_velocity()` の責務を固定
- [ ] **タスク 2.3**: `fix move/surf` での運用確認

### 共通タスク

- [ ] **デバッグフラグの追加**: `DEBUG_SPARTA_SWEEP`, `DEBUG_SPARTA_ROTATE`
- [ ] **コンパイル確認**: エラー・警告なし
- [ ] **単体テスト**: 各関数の動作確認
- [ ] **統合テスト**: 実際のシミュレーションでの動作確認

---

## デバッグとテスト

### デバッグフラグ

```cpp
// src/move_surf.cpp の先頭に追加
#define DEBUG_SPARTA_SWEEP 0   // Phase 1用
#define DEBUG_SPARTA_ROTATE 0  // Phase 2用
```

### テストケース

#### テスト 1: 単純な壁移動
- 2D/3Dの単純な形状（矩形/直方体）を移動
- 壁内粒子が正しく押し出されることを確認

#### テスト 2: 回転壁
- 回転する壁で粒子が正しく反射されることを確認
- `surf_collide` の `translate/rotate` 設定時に期待どおりの反射速度になることを確認

#### テスト 3: INSIDEセル粒子
- INSIDEセル内の粒子が救済または削除されることを確認
- 警告ログが正しく出力されることを確認

### ログ出力例

```
SWEEP: 100 particles detected inside walls
SWEEP: 95 particles swept, 3 rescued, 2 deleted
ROTATE: Reflection handled by surf_collide translate/rotate setting
```

---

## 注意事項

1. **既存コードの保持**: 削除したコードはコメントアウトで保持し、後方互換性を維持
2. **エラーハンドリング**: 射影失敗時は適切にエラー処理を行う
3. **並列処理**: MPI環境での動作を確認（各プロセスで独立に処理）
4. **パフォーマンス**: 近傍探索の効率化を検討（必要に応じて）

---

## 参照ファイル一覧

### 主要ファイル
- `src/move_surf.cpp` - メイン実装
- `src/move_surf.h` - ヘッダー
- `src/fix_move_surf.cpp` - fix move/surf実装

### 関連ファイル
- `src/surf_collide_diffuse.cpp` - Diffuse反射モデル
- `src/surf_collide_specular.cpp` - Specular反射モデル
- `src/surf_collide_cll.cpp` - CLL反射モデル
- `src/grid_surf.cpp` - グリッド-表面関連処理
- `src/geometry.cpp` - 幾何学計算（距離計算等）
- `src/update.cpp` - 粒子更新処理

---

**最終更新**: 2024年  
**ステータス**: 実装前
