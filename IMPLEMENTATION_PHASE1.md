# Phase 1 実装計画：粒子押し出し機能

## 現在の状況

### 既存コードの確認
- `fix_move_surf.cpp`: 壁移動後の処理を実装
- `move_surf.cpp`: 表面移動と粒子削除を実装
- `move_surf::remove_particles()`: 現在は壁内の粒子を**削除**している

### 実装方針
Phase 1では、削除の前に「検知→ログ出力→可視化」を追加し、その後「押し出し処理」を実装する。

## 月1: 検知ロジックの実装

### タスク1.1: 粒子検知機能の追加

**目標**: 壁移動後に「壁内」にある粒子を検知し、ログ出力する機能を実装する。

**実装箇所**: `move_surf.cpp` の `remove_particles()` 関数

**実装内容**:
1. 現在の削除処理の前に、検知ロジックを追加
2. 検知された粒子の情報をログ出力
3. 可視化用のデータ出力（オプション）

**コード構造**:
```cpp
bigint MoveSurf::remove_particles()
{
  // 既存の検知ロジックを確認
  // 削除の前に、検知された粒子の情報を記録
  
  // TODO: Phase 1 月1
  // 1. 壁内粒子の検知（既存ロジックを流用）
  // 2. 検知された粒子の情報をログ出力
  // 3. 可視化用データの出力（dumpファイル等）
  
  // 既存の削除処理（一旦コメントアウトして検証）
  // return ndeleted;
}
```

### タスク1.2: 検証用テストケース

**1次元ピストン問題**:
- シンプルな1次元移動壁
- 粒子が消失する箇所を特定
- ログ出力で消失タイミングを確認

## 月2: 押し出しロジックの実装

### タスク2.1: 基本押し出しアルゴリズム

**実装内容**:
1. **Cell Selection**: 粒子が存在するセルの表面要素を取得
2. **Velocity Filter**: $V_n = \vec{v}_{wall} \cdot \vec{n} > 0$ の面を候補とする
3. **Projection Check**: 射影パラメータ $0 \le t \le 1$ で判定
4. **Re-positioning**: 最短距離の面に安全マージン $\epsilon$ を加えて再配置
5. **Velocity Update**: CLLモデル（または他の表面衝突モデル）を使用して速度を更新
   - `surf->sc[isc]->collide()` を呼び出して反射速度を計算
   - 移動壁を考慮した速度が自動的に計算される

**新しい関数**:
```cpp
// move_surf.h に追加
void push_out_particles();

// move_surf.cpp に実装
void MoveSurf::push_out_particles()
{
  // Step 1: Cell Selection
  // Step 2: Velocity Filter (V_n > 0)
  // Step 3: Projection Check (0 <= t <= 1)
  // Step 4: Re-positioning with epsilon margin
  // Step 5: Velocity Update using CLL model (or other surface collision model)
  //   - surf->sc[isc]->collide() to compute reflection velocity
  //   - Moving wall velocity is automatically considered
}
```

### タスク2.2: 安全マージンの実装

**パラメータ**:
- `epsilon`: 安全マージン（推奨: $10^{-9}$ m）
- 入力ファイルで設定可能にする

## 月3: 並列化対応

### タスク3.1: MPI通信境界の処理

**課題**:
- Ghost cellを跨ぐ粒子移動
- プロセス間での粒子情報の同期

**実装内容**:
- 粒子移動後の通信処理
- メモリリークのチェック

## 実装チェックリスト

### 月1
- [ ] 粒子検知ロジックの実装
- [ ] ログ出力機能
- [ ] 可視化データ出力
- [ ] 1次元ピストンテストケースの作成

### 月2
- [ ] Cell Selection の実装
- [ ] Velocity Filter の実装
- [ ] Projection Check の実装
- [ ] Re-positioning の実装
- [ ] 安全マージンの実装
- [ ] Velocity Update の実装（CLLモデルまたは他の表面衝突モデルを使用）

### 月3
- [ ] MPI通信境界の処理
- [ ] メモリリークチェック
- [ ] ベンチマークテスト（BM-1）

## 関連ファイル

- `src/move_surf.h`: クラス定義
- `src/move_surf.cpp`: 実装
- `src/fix_move_surf.cpp`: Fixクラスのラッパー
- `src/grid.h/cpp`: グリッド関連（set_inout等）

## 次のステップ

1. `remove_particles()` の現在の実装を詳細に確認
2. 検知ロジックの実装開始
3. テストケースの準備




