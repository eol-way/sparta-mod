# 引き継ぎ資料：SPARTA ルーツポンプ移動壁対応（2026年10月15日時点）

この資料は、別の環境で作業を引き継ぐ Claude（または Cursor 等の AI アシスタント）向けです。前の環境での対話履歴は見られない前提で、作業を続けるのに必要な事実・決定事項・注意点をまとめています。

- リポジトリ：`https://github.com/eol-way/sparta-mod.git`
- ブランチ：`phase0-detection`（リモートと同期済みの最新コミットは `164037fc`）
- 設計資料：[PROJECT_PLAN.md](PROJECT_PLAN.md)（Rev2。設計の根拠と全体計画）
- 試験の説明：[examples/moving_wall_validation/README.md](examples/moving_wall_validation/README.md)

---

## 0. 最初に必ず知っておくこと

1. **「テストケース1・2（並進ピストン・二重円筒）の速度分布検証」は実装・実行済みで、全項目合格しています。** 一から作り直さないでください。
   - 作成時点では**未コミット**でした（§4.4 のファイル一覧）。
   - 新しい環境で `git status` を見て、まだ存在するか、コミット済みかを最初に確認してください（§10）。
2. **Phase 番号が2系統あります。**

   | 呼び方 | 内容 |
   |---|---|
   | MVP の作業指示とコード内のタグ | Phase 1 ＝ めり込み補正（幾何学的退避と遅延削除）、Phase 2 ＝ ガリレイ変換による壁速度の伝達 |
   | PROJECT_PLAN.md（Rev2） | Phase 1 ＝ 剛体運動学と壁静止系反射、Phase 2 ＝ めり込み補正 |

   対応表は PROJECT_PLAN.md の §0 にあります。
3. **この作業環境（WSL）には GitHub の認証がありません。** push はユーザーが手動で行っています（VS Code 等）。commit はユーザーの指示があるときだけ行い、メッセージもユーザー指定に従ってきました。

---

## 1. プロジェクトの目的

- オープンソースの DSMC コード SPARTA（Sandia）を改造し、ルーツ真空ポンプの排気を解析できるようにする。
- 対象：2本のローターが逆回転して噛み合い、クリアランスは数 µm。
- 想定スケール：クリアランスの Kn は 1〜1000、翼先端速度は最確速度の 0.06〜0.15 倍、1回転あたり 10⁶〜10⁷ ステップ。
- 上流（未改造）の SPARTA は、移動壁に対して次の致命的な欠陥を持っていた（詳細は PROJECT_PLAN.md §2 の P1〜P14）。
  - **全粒子削除**：`fix move/surf` は、移動した要素を含むセルの粒子を全部消す（`MoveSurf::remove_particles()`）。ローター表面が付着確率 O(1) の吸着面のように振る舞い、試験では 1000 ステップで粒子の 92% が消えた。
  - **壁速度が反射に入らない**：`surf_collide diffuse/cll` の `translate`/`rotate` キーワードは接線成分しか加えず（法線成分を捨てる）、鏡面反射側は実験室系のまま。翼が気体を押す p-dV 仕事が伝わらない。

---

## 2. MVP の成果（すべて `phase0-detection` に push 済み）

| コミット | 内容 |
|---|---|
| `dc4cc6ab` | めり込み補正（`engulf project`）、退避時の反射を壁静止系で実行 |
| `46113cb6` | `Update::move()` の通常の衝突にもガリレイ変換を適用、`wallvel yes\|no`、回帰テスト一式 |
| `be5890b8` | 1ステップ目から正しい壁速度（Nlarge を受け渡し）、`translate/rotate`・`piston` との併用をハードエラー化 |
| `58c633e2` | GIF 自動生成（`generate_all_gifs.sh`、`make_gifs.py`）、GIF 5本 |
| `f8bacab3` | README に GIF を埋め込み |
| `8e8f2d47` | PROJECT_PLAN.md の進捗更新（§0 進捗サマリー） |
| `164037fc` | AMR 格子の回帰テスト（静的細分化、`fix adapt`、3D） |

主な検証結果（GCC + OpenMPI、4 プロセス）：

- **質量保存**：回転楕円（2D）・楕円体（3D）・逆回転2ロータ・`move_surf` コマンドで、粒子数が完全に保存された。全ダンプで壁内の粒子は 0 個。
  - 従来方式では 100000 → 8197 個、3D では 32000 → 4815 個に減る。
- **後方互換**：`engulf delete`（既定）は、未改造の SPARTA と統計出力がビット単位で一致する（`ref/ellipse2d_delete.np4.stats`）。
- **法線方向の運動量伝達**：無衝突気体に 100 m/s で進入する鏡面ピストンの面圧が p/p₀ = 1.4373 で、解析解 1.4377 との差は −0.03%。
  - 解析解：p/p₀ = (1+2s²)(1+erf s) + (2s/√π)e^(−s²)、s = U/c_mp
  - 比較：`wallvel no` では 1.189、従来方式では 0.117。
- **全衝突での運動量伝達**：拡散壁の正方形容器を回すと、気体は剛体回転 u_θ = ωr に到達した（回転率の比 1.0020）。`wallvel no` では比 0.098 にとどまる。
- **起動直後の壁速度**：nevery 20 で最初の 40 ステップを平均したピストン面圧が 1.4386。Nlarge の受け渡しを入れる前の版では 1.3448 だった。
- **めり込み発生率**：掃過体積から求めた理論値 n·(u_w·n)⁺·A と一致した（43.1 個/ステップに対し理論値 42.5）。デバッグビルドで手動確認したもので、回帰テストには含まない。
- **AMR**：静的細分化（最大4段）、`fix adapt`（`surf` と `particle` の両基準）＋`fix balance`、3D の細分化格子で、粒子数の保存と `fix grid/check ... outside yes` の異常件数 0 を確認した。
- **プロセス数**：1・4・8 プロセスのすべてで合格。

---

## 3. 実装した中核技術

### 3.1 前提：時間進行の構造

- 粒子が移動する間、壁は凍結している（`Update::move()`）。
- 壁は `fix move/surf` の `end_of_step()` で、nevery ステップごとに瞬間移動する（ゼロ次ホールド）。
- 移動のたびに格子の全域再カット（`surf2grid`、`set_inout` 等）を行い、その後 `remove_particles()` を呼ぶ。
- このジャンプで壁が掃いた体積の中にいた粒子が「めり込み粒子」になる。

### 3.2 めり込み補正：最近傍射影による退避と遅延削除（`engulf project`）

`MoveSurf::remove_particles()` のセルループの先頭で、`engulf project` なら `push_out_cell()` に分岐する。従来の削除処理は `engulf delete` 用にそのまま残してある。

1. **検出**（`push_out_cell()`）：対象は INSIDE セルの全粒子と、移動した非透過要素を含むセルの粒子。後者は、`Grid::point_outside_surfs()` で得た流体中の参照点までのレイと表面の交差回数の偶奇で、壁内かどうかを判定する（`particle_inside()`。透過面は数えない）。
2. **最近傍の壁面**（`nearest_surf()`）：
   - まず所属セルの表面要素を探す。最短距離が「粒子からセル境界までの距離」以下なら確定する（ほかの要素はセルと重ならないため）。
   - 確定できなければ全要素を探索する（外接箱で枝刈り）。
   - 三角形上の最近点は Ericson の Voronoi 領域法で求め、すべての除算にガードを入れてある。
3. **退避先**（`push_out_one()`）：最近点 xc から `EPSPUSH`（1e-6）×セル寸法だけ流体側へ置く。
   - 方向はまず「粒子 → 最近点」の向き（角部でも有効）、だめなら要素の法線。
   - めり込みが深すぎる粒子（1回の移動量 dmove ＋ セル2個分を超える）は退避しない。
4. **所属セル**（`locate_cell()`）：自セル内でなければ、ルートから `Grid::id_find_child(0,0,boxlo,boxhi,x)` でたどる。
   - 自ランクのセルなら、`point_in_flow()` で流体中かを確認して確定する。
   - ゴーストセル（他ランク）なら移送リストに積み、セルループの後に `Comm::migrate_particles()` で送る。受信側（`receive_pushed()`）で流体中かを確認し、分割セルならサブセルに割り当てる。
5. **反射**（fix の場合のみ）：`reflect_wall_frame()` で壁静止系の反射を行う（§3.3）。`move_surf` コマンドは瞬間移動なので、位置だけ戻して速度は変えない。
6. **削除**：救済できない粒子は `icell = -1` の印付けだけを行う。`remove_particles()` の末尾で `compress_rebalance()` を1回だけ呼んでまとめて消す。
   - 理由別に計上する：no surf / too deep / no cell / blocked / destroyed。
   - 初回だけ警告を出す。fix の出力値（`f_ID`）は削除数になる。

### 3.3 壁静止系反射：ガリレイ変換ラッパー（`wallvel yes`）

- **通常の衝突**：`Update::move()` で、移動表面との衝突ごとに次を行う（目印のコメントは `[ADDED] Phase 2: 移動壁との衝突はガリレイ変換で挟む`）。
  - 壁速度 u_w（衝突点での全ベクトル）を `MoveSurf::wall_velocity_any(sparta, isurf, x, uw)` で得る。
  - `(v − u_w)·n ≥ 0` の場合（`wallframe = 2`）：壁静止系で壁に近づいていない粒子は、後退した壁の古い位置に当たっただけなので、衝突させずにそのまま通過させる。衝突数・表面集計には入れない。
  - それ以外（`wallframe = 1`）：v −= u_w → 既存の `collide()` をそのまま呼ぶ → 元の粒子と、反応で生成された粒子の両方に u_w を加える。
  - 表面集計用の入射状態 `iorig` は実験室系のまま保つ（力・トルクは正しく集計される）。
- **退避時の反射**：`reflect_wall_frame()` で同じ変換を行う。
  - 表面反応は使わない（isr = −1）。粒子が生成されたらエラーにする。
  - `piston` と透過モデルは対象外。
  - この反射も `compute surf` の集計に入れる。
- **有効判定**：`Update::setup()` が run 開始時に `MoveSurf::wallvel_setup()` を1回呼び、`Update::wallvelflag` を設定する。
  - 無効なら衝突ごとの追加処理は整数判定1回だけで、結果は従来と同一。
  - `wallvel yes` の表面に、`piston` か、`translate`/`rotate` 付きの `diffuse`/`cll` が割り当てられていれば、`error->all()` で停止する。判定には `SurfCollide::has_added_velocity()`（基底では false、diffuse/cll が `trflag` を返す）を使う。
- **移動グループの一覧**：`move_surf.cpp` 内の static 変数 `moving_groups`。`wallvel yes` の `fix move/surf` が `process_args()` で登録し、デストラクタで外す。

### 3.4 壁速度の算出

- `FixMoveSurf` のコンストラクタが `movesurf->set_nlarge(nlarge)` を、`FixMoveSurf::init()`（run ごと）が `movesurf->init()` を呼ぶ。
- 回転は ω = θ/(Nlarge·dt)（`omega` × `rvec` で回転軸方向）、並進は速度 = Δ/(Nlarge·dt)。run の1ステップ目から有効。
- 表面が動くたび（`update_kinematics()`）にも同じ式で再計算するので、timestep の変更に追従する。

### 3.5 入力キーワードと制約

| キーワード | 値 | 既定 |
|---|---|---|
| `engulf` | `delete`（従来）/ `project` | `delete` |
| `wallvel` | `yes` / `no` | `project` なら yes、`delete` なら no。yes には project が必須 |

入力例：

```
fix 10 move/surf rotor 1 2000 rotate 360 0 0 1 5 5 0 engulf project
```

`engulf project` では次をエラーにしている：Kokkos、implicit surf。

### 3.6 変更したファイル

| ファイル | 内容 |
|---|---|
| `src/move_surf.cpp/.h` | めり込み補正、壁速度の算出と参照、併用チェック（上記の関数群） |
| `src/update.cpp/.h` | 衝突処理へのガリレイ変換、`wallvelflag` |
| `src/fix_move_surf.cpp` | Nlarge の受け渡し、run ごとの初期化 |
| `src/surf_collide.h`、`surf_collide_diffuse.h`、`surf_collide_cll.h` | `has_added_velocity()` |

コードの約束事：

- 変更箇所には `// [ADDED] Phase X: 説明` / `[MODIFIED]` / `[REMOVED]`（削除したコードはコメントで残す）を付ける。番号は MVP の作業指示の Phase。
- デバッグ出力は `DEBUG_SPARTA_ENGULF` / `DEBUG_SPARTA_WALLFRAME`。`#ifndef` 形式なので、`-D...=1` を付けてコンパイルすると有効になる（ステップごとの件数と ω を表示）。

---

## 4. 直近で合意した検証方針：速度分布の検証（気体衝突なし）

### 4.1 目的

移動壁での拡散反射が、壁静止系で正しく行われていることを、反射直後の粒子の速度分布そのもので示す。気体分子どうしの衝突は切り（`collide` コマンドなし）、壁と粒子の相互作用だけを測る。

### 4.2 設計（実装済み）

| | テストケース1：並進（`in.test_translation`） | テストケース2：回転（`in.test_rotation`） |
|---|---|---|
| 次元・形状 | 3D の管路 x∈[0,3]、y,z∈[0,0.2]。ピストンは x=0.2 の三角形2枚（`data.piston3d`、法線 +x） | 2D。内円 R=0.5（`data.ring_inner`、72 本、法線は外向き）、外円 R=1.0（`data.ring_outer`、72 本、法線は内向き） |
| 格子 | 60×4×4、境界 `o r r` | 43×43（セル幅 0.0512）、境界 `o o p` |
| 壁の動き | `trans 1.0 0 0` を 1000 ステップ、dt 1e-5 → U = 100 m/s | 内円を 1600 ステップで1回転 → ω = 392.7 rad/s、Rω = 196.3 m/s |
| 反射 | `diffuse 300.0 1.0`（対照は `vanish`） | 内円は同左、外円は `vanish` |
| 入射粒子 | `fix emit/face beam xhi`（vstream −1000、30 K のビーム） | `fix emit/surf beam outer normal yes`（中心向き 1000 m/s、30 K） |
| 出力 | step 400〜1000、10 ステップごと、x < 1.3 の粒子（id x vx vy vz） | step 200〜600、10 ステップごと、r < 0.6 の粒子（id x y z vx vy vz） |
| 抽出 | x_p < x < x_p+0.05 かつ w = vx − U > 0（x_p = 0.2 + step/1000） | 0.5 < r < 0.55 かつ v_r > 0 |

設計上の要点（いずれも検討のうえ採用したもの）：

- **入射をビームにした理由**：入射粒子がすべて壁に向かうので、「壁静止系で離れていく粒子」は反射粒子だけになる。壁を `vanish` にした対照実行で、この判定に紛れ込む粒子が 0 個であることを実測している。拡散反射の出射分布は入射分布に依存しないので、ビームにしても検証の一般性は失われない。
- **密度重み付けから流束重み付けへの換算**：層の中の粒子は遅いものほど長く層にとどまるので、そのままでは密度重み付きの分布になる。各粒子を層内の滞在時間の逆数で重み付けし、壁から出ていく流束の分布（散乱核そのもの）に換算して比べる。
  - 並進：滞在時間は 0.05/w_n
  - 回転：直線軌道で r = R から r = 0.55 に達するまでの時間
- **円筒の曲率補正**：層の中では周方向速度が角運動量保存で R/r 倍に下がる。各粒子を直線軌道で r = R の放出点までさかのぼり、その点の法線・接線方向で成分を求める。
- **回転を 2D にした理由**：3D の円筒は、z 面に接したまま回すと上流の `rotate_3d()` の丸め誤差で端の頂点が箱の外に出て止まる（§6）。

### 4.3 理論と結果

理論（c = √(2kT_w/m) = 596.9 m/s、N 原子、T_w = 300 K）：

- 法線成分：(2w/c²)exp(−w²/c²) のレイリー分布。平均 √π·c/2 = 529.0、二乗平均平方根 c。
- 接線成分：ガウス分布。平均は壁速度の接線成分（並進 0、回転 196.3）、標準偏差 c/√2 = 422.1。
- 合否の基準：平均・標準偏差が「標準誤差の4倍＋c の 1%」以内。吸収壁の対照で、反射粒子と誤認される粒子が 0.1% 以下。

| 量 | 並進 `wallvel yes` | 回転 `wallvel yes` | 理論 | `wallvel no` |
|---|---|---|---|---|
| 法線成分の平均 | 530.5 | 528.3 | 529.0 | 並進で 445.3 |
| 接線成分の平均 | 1.3 | 196.9 | 0 / 196.3 | 回転で 0.5 |
| 接線成分の標準偏差 | 421.8 | 421.8 | 422.1 | 一致 |

- 14 項目すべて合格。反射粒子のサンプル数は並進で約 40 万、回転で約 32 万。
- 吸収壁の対照で、離れていく粒子は両ケースとも 0 個。
- 図：`examples/moving_wall_validation/plots/velocity_translation.png`、`velocity_rotation.png`

### 4.4 該当ファイル（作成時点で未コミット）

`examples/moving_wall_validation/` の次のファイルです。

- 新規：`in.test_translation`、`in.test_rotation`、`make_surfs.py`（形状ファイルの生成）、`data.piston3d`、`data.ring_inner`、`data.ring_outer`、`verify_distribution.py`、`plots/velocity_*.png`
- 変更：`run_all.sh`（上記6ケースと `verify_distribution.py` の呼び出しを追加）、`README.md`（「6. 拡散反射の速度分布」の節）

`verify_distribution.py` の判定は純 Python で動きます。matplotlib があれば図も出力します。

使い方：`python3 verify_distribution.py OUTDIR [PLOTDIR]`

---

## 5. 回帰テストと可視化の使い方

```
cd src && make -j 20 mpi                          # src/spa_mpi を生成（20並列で数秒）
cd ../examples/moving_wall_validation
./run_all.sh [実行ファイル] [プロセス数]           # 既定 ../../src/spa_mpi、4
./generate_all_gifs.sh                            # GIF を gifs/ に再生成（約 30 秒）
```

- `run_all.sh` は `out/`（git 管理外）で実行し、`check.py` と `verify_distribution.py` で合否を判定する。
  - 4 プロセスで約 1 分。期待値は既存 39 項目＋速度分布 14 項目の合格（1・8 プロセスでは参照比較が SKIP になる）。
  - 未改造版との一致比較の参照は 4 プロセス専用。
- 試験は全部で、回転楕円（2D）、楕円体（3D）、2ロータ、`move_surf` コマンド、ピストン（較正・壁速度あり／なし・従来・起動直後）、回転容器、AMR（6ケース）、併用禁止の確認（3件は停止を期待）、速度分布（6ケース）。
- GIF・README の図は git で管理する（ユーザーの決定）。

---

## 6. 既知の問題・注意点（作業中に判明したもの）

- **形状の退化配置（上流の問題）**：面や頂点が格子線・格子点に正確に乗ると、`Grid::set_inout()` が「Cell type mis-match」で停止する。従来モードでも同じ。試験形状は格子からずらしている（例：`data.square` の一辺を 3.9、回転試験の格子を 43×43）。
- **3D の回転（上流の問題）**：箱の面に接する開いた 3D 表面を `rotate` すると、回転行列の丸め誤差で頂点が箱の外に出て止まる。端面を持つ閉じた形状なら当たらない。
- **dump の出力精度**：既定の 6 桁では、壁から セル幅×10⁻⁶ だけ離れた粒子が丸めで「壁内」に見える。内外判定に使う dump は `dump_modify ... format float %20.15g` にする。
- **SPARTA の文法上の注意**
  - `dump_modify sort` は無い。
  - `compute surf` は値が1つでも per-surf 配列なので `c_ID[1]` と書く。
  - `stats` の出力間隔は `fix ave/*` の Nfreq と合わせる。
  - `fix adapt` には `iterate` を書けない。
  - `fix adapt ... surf` はセル寸法でしか粗大化しないので、壁が通過した跡に細かいセルが残る。
- **性能**：最寄り壁面を確定できない場合の全要素探索は O(N_surf)。3D の本番規模では遅い（PROJECT_PLAN.md Phase 4 の課題）。
- **未対応**：Kokkos、implicit surf、退避時の表面反応。計画書 Phase 2 の遡及衝突（鏡像法）・統計的再放出、Phase 3（解析 SDF）、Phase 4（連続衝突判定）は未着手。

---

## 7. 作業環境（前の環境）

- WSL2 Ubuntu、リポジトリは `/home/yu26016/SPARTA/GIT/sparta-mod`。
- ビルド：`src/` で `make -j 20 mpi`（`MAKE/Makefile.mpi`、OpenMPI の `mpicxx`）。`Obj_mpi/` と実行ファイルは git 管理外。
- 実行：`mpirun --oversubscribe -np N src/spa_mpi -in ...`
- Python 3.12。システムには numpy・matplotlib が無い。`examples/moving_wall_validation/.venv`（git 管理外）に numpy・matplotlib・Pillow を入れて使っている。無い環境では次で作る。

  ```
  python3 -m venv .venv && .venv/bin/pip install matplotlib pillow
  ```

---

## 8. ユーザーとの作業の進め方

- 日本語でやり取りする。設計の議論と実装を分けて進める方針。
- 他の AI（Gemini 等）が書いたプロンプトを貼り付けて依頼することが多い。技術的な前提は鵜呑みにせず確認し、変えた点は理由とともに報告する。
  - 例：直近の依頼では、反射直後の分布を理論のマクスウェル分布と直接比べる想定だったが、密度重み付けと流束重み付けの違いを補正した。
- 依頼文で編集してよいファイルが指定されることがある。範囲を守り、範囲内で目的を達成できない場合はその理由を明示する。
- commit・push の流儀は §0 の 3 のとおり（commit は指示があるときだけ、push はユーザー）。
- 結果は数値と根拠（実行ログ）で報告する。ビルドが通っただけで完了とせず、動作を検証する。

---

## 9. 決まっていないこと

- PROJECT_PLAN.md §12 の未決事項はすべて未確定。
  - 代表クリアランス h（試験では 5 µm を仮定）
  - ローター翼形の定義元（解析式か CAD か）
  - 2D 先行か 3D か
  - 壁の反射モデルと適応係数
  - GPU が必須か
- MVP 後の次の開発テーマ（下記の候補のどれから進めるか）は、ユーザーがまだ決めていない。

---

## 10. Next Action（新しい Claude が最初に行うこと）

1. **状態の確認**
   - `git status` と `git log --oneline -3` を見る。
   - §4.4 のファイルが作業ツリーにあり、未コミットなら、ユーザーにコミットしてよいか確認する。メッセージはユーザーの指定に従う。
   - ファイルが存在しない（前の環境でコミット・push されずに移行した）場合は、§4.2〜4.3 の仕様で再作成する。再作成前に、前の環境に残っていないかユーザーに確認する。
2. **ビルドと回帰テスト**：§5 の手順で `make` → `./run_all.sh` を実行し、全項目が合格することを確かめる（新しい環境の健全性確認）。
3. **次のタスクをユーザーに確認する**。依頼文に「テストケース1・2の Python 検証スクリプト作成」があれば、すでに完了していることを伝える。候補（いずれも未決定）：
   - 速度分布の検証を、鏡面反射・Maxwell 型（acc < 1）・CLL に広げる（同じ枠組みで理論分布を書ける）。
   - PROJECT_PLAN.md の Phase 3（解析 SDF）・Phase 4（連続衝突判定と高速化）。
   - Kokkos・implicit surf・表面反応への対応。
   - ルーツポンプの実形状（2本の翼形とケーシング）での試計算。
