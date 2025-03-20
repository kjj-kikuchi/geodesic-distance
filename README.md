# 測地距離の計算

<img src = "kitten_geodesic1.png" width = 50%><img src = "kitten_geodesic2.png" width = 50%>  

三角形メッシュの測地距離の計算手法 Geodesics in heat: A new approach to computing distance based on heat flow [Crane et al. 2013]の実装．  

* 使用言語：C++
* 使用ライブラリ：Eigen


### サンプル出力の情報
ソースとなる頂点を頂点7000とした結果．
- 入力メッシュ : `294_kitten_uniform.obj` （頂点数 134,448 面数 268,896）
- 出力：`kitten_uniform_geodesic.vtk`
- 実行時間 : 2.72206 秒
