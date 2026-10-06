# moving_wall_validation — 移動壁の回帰テスト

`fix move/surf` / `move_surf` に追加した2つの機能を検証する回帰テストです。

| キーワード | 値 | 内容 |
|---|---|---|
| `engulf` | `delete`（既定）/ `project` | `delete` = 従来どおり、移動要素を含むセルの全粒子を削除。`project` = 壁内に入った粒子だけを最寄りの壁面の外側へ退避（他ランクのセルへは粒子移送）、救済できない粒子のみ計上して削除 |
| `wallvel` | `yes` / `no` | `yes` = 移動表面とのすべての衝突（通常の粒子移動中と退避時）を壁静止系で行う：`w = v − u_w` → 既存の `collide()` → `v' = w' + u_w`。壁速度は ω = θ/(Nlarge·dt)、並進速度 = Δ/(Nlarge·dt) で、run の1ステップ目から有効。既定は `engulf project` なら `yes`、`delete` なら `no`。`yes` は `engulf project` が必須 |

入力例：

```
fix 10 move/surf rotor 1 2000 rotate 360 0 0 1 5 5 0 engulf project
```

## 実行方法

```
./run_all.sh [SPARTA実行ファイル] [プロセス数]
# 既定: ../../src/spa_mpi, 4 プロセス
# 起動コマンドは環境変数 MPIRUN で変更可（既定: mpirun --oversubscribe -np N）
```

出力はすべて `out/` に書かれ（git 管理外）、最後に `check.py` が各項目を `PASS` / `FAIL` で判定します。1つでも `FAIL` があるか、実行が異常終了すると終了コード 1 を返します。所要時間は 4 プロセスで約 20 秒です。

## 試験一覧

| 入力 | 検証内容 | 合格条件 | 実測（GCC + OpenMPI, 4 プロセス） |
|---|---|---|---|
| `in.ellipse2d` | 閉じた箱の中で回転する楕円（2D）。`delete` と `project` を比較 | `project`：粒子数が不変、どの時刻でも壁内の粒子 0 個。`delete`：未改造 SPARTA と統計出力が一致（`ref/`、4 プロセス時のみ） | 100000 → 100000、壁内 0。`delete` は 100000 → 8197 で参照と完全一致 |
| `in.ellipsoid3d` | 回転する楕円体（3D、三角形要素） | 粒子数が不変、壁内 0 個 | 32000 → 32000、壁内 0（`delete` は → 4815） |
| `in.two_rotors` | 逆回転する2本の楕円（ルーツ型、fix 2本） | 粒子数が不変、両ロータの内部 0 個 | 100000 → 100000、内部 0 |
| `in.move_surf_cmd` | `move_surf` コマンドによる瞬間移動（セル4個分の平行移動） | 削除 0、移動後の内部 0 個 | 退避 6688、削除 0（`delete` は 8596 個削除） |
| `in.piston` | 鏡面反射のピストンが U = 100 m/s で無衝突気体に進入。法線方向の運動量伝達を検証 | 静止壁の較正 p/p₀ = 1 ± 3%。`wallvel yes` の面圧が自由分子流の解析解 p/p₀ = (1+2s²)(1+erf s) + (2s/√π)e^(−s²), s = U/c_mp の ±3% 以内 | 較正 1.0039、`wallvel yes` 1.4373（解析解 1.4377、差 −0.03%）。参考：`wallvel no` 1.1887、`delete` 0.1166 |
| `in.piston`（初動） | `nevery 20` で最初の 40 ステップだけを平均。最初の表面移動より前から壁速度が正しく入るかを検証 | 解析解の ±3% 以内 | 1.4386（+0.06%）。壁速度を初回移動時に測っていた旧版では 1.3448（−6.5%） |
| `in.rotating_box` | 静止した気体を入れた正方形容器が ω = 62.83 rad/s で回転（拡散反射、壁温＝気体温度）。接線方向を含む全衝突での運動量伝達を検証 | 粒子数が不変。平衡解である剛体回転 u_θ = ωr に対し、気体の回転率（r < 1.75 のセルで最小二乗）が ω の ±3% 以内 | 回転率 62.96 1/s（比 1.0020）。参考：`wallvel no` では比 0.0975 |
| `in.guard_rotate` / `in.guard_cll_translate` / `in.guard_piston` | 誤った組み合わせの拒否：`wallvel yes` の表面に、自分で壁速度を加える反射モデル（`diffuse` の `rotate`、`cll` の `translate`、`piston`）を割り当てる | `error->all()` で停止し、所定のエラーメッセージが出ること。`wallvel no` との併用は従来どおり動くこと | 3件とも停止、`wallvel no` は実行可 |

`INFO` 行は比較用の参考値で、合否には使いません。従来方式や `wallvel no` では検証量が大きく外れることを示しています。

## 注意点

- **参照データ**：`ref/ellipse2d_delete.np4.stats` は、未改造 SPARTA（commit 62cc4756）を GCC + OpenMPI で 4 プロセス実行して作ったものです。プロセス数が異なると比較は `SKIP`、別のコンパイラや MPI で浮動小数点の差が出た場合は `WARN` になります。
- **`data.square`**：辺を格子線（セル幅 0.25）から外すため、一辺を 3.9 にしています。面がセルの面とちょうど重なると、SPARTA 本体の `Grid::set_inout()` が「Cell type mis-match」で停止します（従来の `engulf delete` でも同じです）。
- **壁速度の初期化**：`fix move/surf` が `Nlarge` を渡し、run ごとの `init()` でその時点の timestep から壁速度を計算します。表面の移動のたびにも同じ式で再計算するので、run の途中で timestep が変わっても追従します。
- **反射モデルの組み合わせ**：`wallvel yes` の表面に `surf_collide piston`、または `translate` / `rotate` キーワード付きの `diffuse` / `cll` を割り当てると、run 開始時にエラーで停止します（壁速度の二重加算を防ぐため。判定は `SurfCollide::has_added_velocity()`）。
- **未対応**：Kokkos と implicit surf では `engulf project` をエラーにしています。
- **詳細な計数**：`src/move_surf.cpp` を `-DDEBUG_SPARTA_ENGULF=1 -DDEBUG_SPARTA_WALLFRAME=1` 付きでコンパイルすると、ステップごとに壁内粒子・退避・壁静止系での反射の件数と、ω・最大移動量が表示されます。

## GIF アニメーション

```
./generate_all_gifs.sh [SPARTA実行ファイル] [プロセス数]
```

入力を `-var anim N`（N ステップごとに粒子のサンプルと壁形状を dump）で実行し、`make_gifs.py` が `gifs/` に GIF を書き出します。作業ファイルは `anim_out/`（git 管理外）に置かれ、所要時間は 4 プロセスで約 30 秒です。`-var anim` を付けない通常実行では dump は出力されず、回帰テストの結果は変わりません。

Python には numpy・matplotlib・Pillow が必要です。`PYTHON` 環境変数、`./.venv/bin/python`、`python3` の順に探し、無ければ次の導入手順を表示して停止します。

```
python3 -m venv .venv && .venv/bin/pip install matplotlib pillow
```

| GIF | 内容 |
|---|---|
| `gifs/ellipse2d.gif` | 回転楕円。左が従来の `engulf delete`（粒子が消えていく）、右が `engulf project`（粒子数を維持） |
| `gifs/two_rotors.gif` | 逆回転する2本のロータ（`engulf project`） |
| `gifs/ellipsoid3d.gif` | 回転する3D楕円体。中央の層 \|z\| < 0.2 の粒子と z = 0 断面。左が `delete`、右が `project` |
| `gifs/piston.gif` | 進入するピストン。上3段が `wallvel yes` / `wallvel no` / `delete` の粒子、下段が密度分布 n/n₀（`wallvel yes` ではピストン前方に圧縮領域ができる） |
| `gifs/rotating_box.gif` | 回転する正方形容器。下段は平均周方向速度 u_θ(r)。`wallvel yes` は剛体回転 ωr に到達、`wallvel no` はほぼ静止のまま |

壁の内側（固体側）にある粒子は赤で大きく描き、各フレームの見出しに件数（inside wall）を表示します。固体側の判定は、2D では法線の向き（SPARTA の法線 = ẑ×(v2−v1) が気体側）を考慮した巻き数、3D では凸多面体の面判定で行います。dump は全精度（%20.15g）で出力するので、壁から 10⁻⁶ セル幅だけ離して置いた粒子が丸めで内側と判定されることはありません。現在の結果では、全 GIF の全フレームで 0 件です。
