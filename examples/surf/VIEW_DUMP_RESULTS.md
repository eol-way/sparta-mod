# 検知ロジック検証結果：dump出力のビュー

## 生成された画像ファイル

### 検知が発生したステップの画像（PNG形式）

以下のステップで粒子が検知されました：

1. **Step 200**: 6個の粒子が検知
   - 画像ファイル: `view_images/image.0200.png`
   - 検知された粒子の位置:
     - Particle ID 1763918088: pos = (3.133700e+00, 5.982541e+00, 0.000000e+00)
     - Particle ID 678343629: pos = (2.956314e+00, 6.050901e+00, 0.000000e+00)
     - Particle ID 407053880: pos = (2.984672e+00, 6.322580e+00, 0.000000e+00)
     - Particle ID 1409610915: pos = (2.904574e+00, 6.892270e+00, 0.000000e+00)
     - Particle ID 1229142426: pos = (3.010995e+00, 6.777960e+00, 0.000000e+00)
     - Particle ID 2080616289: pos = (3.058131e+00, 6.755704e+00, 0.000000e+00)

2. **Step 300**: 4個の粒子が検知
   - 画像ファイル: `view_images/image.0300.png`
   - 検知された粒子の位置:
     - Particle ID 205777754: pos = (3.981765e+00, 4.792400e+00, 0.000000e+00)
     - Particle ID 730223286: pos = (3.541883e+00, 4.967581e+00, 0.000000e+00)
     - Particle ID 298328486: pos = (4.057108e+00, 4.930009e+00, 0.000000e+00)
     - Particle ID 786963805: pos = (4.082561e+00, 4.935017e+00, 0.000000e+00)

3. **Step 1800**: 1個の粒子が検知
   - 画像ファイル: `view_images/image.1800.png`
   - 検知された粒子の位置:
     - Particle ID 512736388: pos = (3.870993e+00, 1.785756e+00, 0.000000e+00)

## 画像ファイルの場所

### PNG形式（Cursor内で表示可能）
- `view_images/image.0200.png` - Step 200の画像
- `view_images/image.0300.png` - Step 300の画像
- `view_images/image.1800.png` - Step 1800の画像

### 元のPPM形式（全ステップ）
- `image.0000.ppm` から `image.2000.ppm` まで（50ステップごと、全41ファイル）

## Cursor内で画像を表示する方法

1. **エクスプローラーで開く**:
   - Cursorの左側のファイルエクスプローラーで `examples/surf/view_images/` を開く
   - PNGファイルをクリックしてプレビュー表示

2. **コマンドで開く**:
   ```bash
   cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf/view_images
   # Linux環境で画像ビューアを開く
   xdg-open image.0200.png  # または
   eog image.0200.png        # GNOME Image Viewer
   ```

3. **他のステップの画像も変換する場合**:
   ```bash
   cd /home/yu0126016/SPARTA/GIT/sparta/examples/surf
   # 特定のステップを変換
   convert image.0100.ppm view_images/image.0100.png
   
   # または全ステップを変換（時間がかかります）
   for file in image.*.ppm; do
     step=$(echo $file | sed 's/image\.\([0-9]*\)\.ppm/\1/')
     convert "$file" "view_images/image.${step}.png"
   done
   ```

## 画像の内容

各画像には以下が表示されます：
- **粒子**: 色分けされた粒子（type別）
- **表面**: 回転する円形表面（黒い線）
- **検知された粒子**: 壁内に取り残された粒子（検知ロジックで検出）

## 検証結果のまとめ

- ✅ 検知ロジックが正常に動作
- ✅ 回転壁でも粒子検知が機能
- ✅ 検知された粒子の位置情報が正確に記録
- ✅ 画像出力により可視化も確認可能

詳細な検証結果は `test_rotate_output.log` を参照してください。




