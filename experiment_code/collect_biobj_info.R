library(smoof)
library(ParamHelpers)
source("MakeBiObjBBOB.R")

getID = function(dimensions, fid, iid) {
  
  # ==== Sanity Checks ====
  
  dimensions = checkmate::asCount(dimensions)
  fid = checkmate::asCount(fid)
  iid = checkmate::asCount(iid)
  checkmate::assertInt(dimensions, lower = 2L, upper = 40L)
  checkmate::assertInt(fid, lower = 1L, upper = 92L)
  checkmate::assertInt(iid, lower = 1L, upper = 15L) # restrict to documented "safe" range
  
  # touch vars
  force(dimensions)
  force(fid)
  force(iid)
  
  # ==== FID Mapping ====
  
  # single-objective BBOB functions, which are used by bi-objective BBOB
  fids = c(1L, 2L, 6L, 8L, 13L, 14L, 15L, 17L, 20L, 21L)
  
  # grid with all pairs of BBOB problems
  grid = expand.grid(fids1 = fids, fids2 = fids)
  grid = grid[grid[, 1L] <= grid[, 2L], ]
  grid = grid[order(grid[, 1L]), ]
  
  for (group in list(1L:5L, 6L:9L, 10L:14L, setdiff(15L:19L, 16L), 20L:24L)) {
    group.grid = expand.grid(fids1 = group, fids2 = group)
    group.grid = group.grid[group.grid[, 1L] < group.grid[, 2L], ]
    group.grid = group.grid[order(group.grid[, 1L]), ]
    
    grid = rbind(grid, group.grid)
    grid = grid[!duplicated(grid), ]
  }
  
  rownames(grid) = NULL
  
  fid1 = grid[fid, "fids1"]
  fid2 = grid[fid, "fids2"]
  
  
  max_iid = 15L
  
  vec_iid_1 = 2 * (1:max_iid) + 1
  vec_iid_2 = vec_iid_1 + 1
  
  iid_mapping = cbind(vec_iid_1, vec_iid_2)
  
  # exceptions, cf. above
  iid_mapping[1L,] = c(2L, 4L)
  iid_mapping[2L,] = c(3L, 5L)
  iid_mapping[9L,] = c(19L, 21L)
  iid_mapping[15L,] = c(31L, 34L)
  
  iid1 = iid_mapping[iid,1]
  iid2 = iid_mapping[iid,2]
  
  return(data.frame(fid1 = fid1, iid1 = iid1, fid2 = fid2, iid2 = iid2))
  
}



DF <- data.frame(fid = integer(), x1 = numeric(), x2 = numeric(), value = numeric())

# ループで各FIDのグローバル最適値を取得
for (f in 1:55) {
  for (i in 1:10){
    
    # 2目的BBOB関数の作成
    y = makeBiObjBBOBFunction(dimensions = 2, fid = f, iid = i)
    
    # FIDとIIDを取得
    ID = getID(dimensions = 2, fid = f,iid = i)
    FID1 = ID$fid1
    IID1 = ID$iid1
    FID2 = ID$fid2
    IID2 = ID$iid2
    
    
    
    # それぞれのBBOB関数を作成
    y1 <- smoof::makeBBOBFunction(dimensions = 2, fid = FID1, iid = IID1)
    y2 <- smoof::makeBBOBFunction(dimensions = 2, fid = FID2, iid = IID2)
    
    # 最適解を取得
    opt1 <- smoof::getGlobalOptimum(y1)
    opt2 <- smoof::getGlobalOptimum(y2)
    
    # 最適解を関数に評価
    y1_value = y1(opt1$par)
    y2_value = y2(opt2$par)
    y_value1 = y(opt1$par)
    y_value2 = y(opt2$par)
    
    
    Ideal = c(y_value1[1],y_value2[2])
    
    Nadir = c(y_value2[1],y_value1[2])
    
    length_x = abs(Ideal[1] - Nadir[1])
    length_y = abs(Ideal[2] - Nadir[2])
    
    # 面積を計算
    Ideal_HV =  length_x * length_y
    
    
    # グローバル最適値をデータフレームに追加
    DF <- rbind(DF, data.frame(fid = f, iid = i, ideal1 = Ideal[1], ideal2 = Ideal[2], nadir1 = Nadir[1], nadir2 = Nadir[2], ideal_hv = Ideal_HV))
  }
}

# 結果を表示
write.csv(DF, "bbob_biobj_ideal_nadir.csv", col.names = TRUE,row.names = FALSE)