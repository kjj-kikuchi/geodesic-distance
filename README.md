# 測地距離の計算

<img src = "kitten_geodesic1.png" width = 50%><img src = "kitten_geodesic2.png" width = 50%>  

三角形メッシュの測地距離の計算手法 Geodesics in heat: A new approach to computing distance based on heat flow [Crane et al. 2013]の実装．
### 概略
1. `make_cotan_laplacian()`でcotan重みラプラシアン`L_cotan`を作る．
1. `make_mass_matrix()`で質量行列Aを作る．
1. `get_time()`でパラメータ`time`の値を得る．
1. ラプラシアンの線形システムを解いて熱流関数uを計算する．

1. `get_gradient()`でuの勾配∇uを計算する．
1. `get_divergence()`でベクトル場Xの発散∇・Xを計算する．

1. ラプラシアンの線形システムを解いて測地距離を計算する．

### サンプル出力の情報
ソースとなる頂点を頂点7000とした結果．
- 入力メッシュ : `294_kitten_uniform.obj` （頂点数 134,448 面数 268,896）
- 出力：`kitten_uniform_geodesic.vtk`
- 実行時間 : 2.72206 秒
