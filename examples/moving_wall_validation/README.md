# moving_wall_validation — 移動壁の回帰テスト

`fix move/surf` / `move_surf` に追加した2つの機能を検証する回帰テストです。

| キーワード | 値 | 内容 |
|---|---|---|
| `engulf` | `delete`（既定）/ `project` | `delete` = 従来どおり、移動要素を含むセルの全粒子を削除。`project` = 壁内に入った粒子だけを最寄りの壁面の外側へ退避（他ランクのセルへは粒子移送）、救済できない粒子のみ計上して削除 |
| `wallvel` | `yes` / `no` | `yes` = 移動表面とのすべての衝突（通常の粒子移動中と退避時）を壁静止系で行う：`w = v − u_w` → 既存の `collide()` → `v' = w' + u_w`。既定は `engulf project` なら `yes`、`delete` なら `no`。`yes` は `engulf project` が必須 |

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

出力はすべて `out/` に書かれ（git 管理外）、最後に `check.py` が各項目を `PASS` / `FAIL` で判定します。1つでも `FAIL` があるか、実行が異常終了すると終了コード 1 を返します。所要時間は 4 プロセスで約 15 秒です。

## 試験一覧

| 入力 | 検証内容 | 合格条件 | 実測（GCC + OpenMPI, 4 プロセス） |
|---|---|---|---|
| `in.ellipse2d` | 閉じた箱の中で回転する楕円（2D）。`delete` と `project` を比較 | `project`：粒子数が不変、どの時刻でも壁内の粒子 0 個。`delete`：未改造 SPARTA と統計出力が一致（`ref/`、4 プロセス時のみ） | 100000 → 100000、壁内 0。`delete` は 100000 → 8197 で参照と完全一致 |
| `in.ellipsoid3d` | 回転する楕円体（3D、三角形要素） | 粒子数が不変、壁内 0 個 | 32000 → 32000、壁内 0（`delete` は → 4815） |
| `in.two_rotors` | 逆回転する2本の楕円（ルーツ型、fix 2本） | 粒子数が不変、両ロータの内部 0 個 | 100000 → 100000、内部 0 |
| `in.move_surf_cmd` | `move_surf` コマンドによる瞬間移動（セル4個分の平行移動） | 削除 0、移動後の内部 0 個 | 退避 6688、削除 0（`delete` は 8596 個削除） |
| `in.piston` | 鏡面反射のピストンが U = 100 m/s で無衝突気体に進入。法線方向の運動量伝達を検証 | 静止壁の較正 p/p₀ = 1 ± 3%。`wallvel yes` の面圧が自由分子流の解析解 p/p₀ = (1+2s²)(1+erf s) + (2s/√π)e^(−s²), s = U/c_mp の ±3% 以内 | 較正 1.0039、`wallvel yes` 1.4368（解析解 1.4377、差 −0.06%）。参考：`wallvel no` 1.1887、`delete` 0.1166 |
| `in.rotating_box` | 静止した気体を入れた正方形容器が ω = 62.83 rad/s で回転（拡散反射、壁温＝気体温度）。接線方向を含む全衝突での運動量伝達を検証 | 粒子数が不変。平衡解である剛体回転 u_θ = ωr に対し、気体の回転率（r < 1.75 のセルで最小二乗）が ω の ±3% 以内 | 回転率 62.88 1/s（比 1.0008）。参考：`wallvel no` では比 0.0975 |

`INFO` 行は比較用の参考値で、合否には使いません。従来方式や `wallvel no` では検証量が大きく外れることを示しています。

## 注意点

- **参照データ**：`ref/ellipse2d_delete.np4.stats` は、未改造 SPARTA（commit 62cc4756）を GCC + OpenMPI で 4 プロセス実行して作ったものです。プロセス数が異なると比較は `SKIP`、別のコンパイラや MPI で浮動小数点の差が出た場合は `WARN` になります。
- **`data.square`**：辺を格子線（セル幅 0.25）から外すため、一辺を 3.9 にしています。面がセルの面とちょうど重なると、SPARTA 本体の `Grid::set_inout()` が「Cell type mis-match」で停止します（従来の `engulf delete` でも同じです）。
- **壁速度の立ち上がり**：壁速度は表面が初めて動いた時点（`nevery` ステップ目）の移動量から求めます（ω = θ·Δfraction / (Δstep·dt)）。それより前の `nevery` ステップは壁速度 0 として扱います。
- **反射モデルの組み合わせ**：`wallvel yes` の表面に `surf_collide piston` を割り当てるとエラーになります。`diffuse` / `cll` の `translate` / `rotate` キーワードは検出できず、併用すると壁速度が二重に加算されるので使わないでください。
- **未対応**：Kokkos と implicit surf では `engulf project` をエラーにしています。
- **詳細な計数**：`src/move_surf.cpp` を `-DDEBUG_SPARTA_ENGULF=1 -DDEBUG_SPARTA_WALLFRAME=1` 付きでコンパイルすると、ステップごとに壁内粒子・退避・壁静止系での反射の件数と、ω・最大移動量が表示されます。
