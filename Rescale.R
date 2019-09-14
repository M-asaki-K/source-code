#Re-scale of optimized x
A <- as.matrix(t(apply(x[ , ], 2, sd)))
B <- as.matrix(t(apply(x[ , ], 2, mean)))
#preprocessed.xの部分は規格化された任意のデータセットに変更可能（例えば逆解析で出た最適条件とか）
xrev <- t(apply(preprocessed.x, 1, function(preprocessed.x){(preprocessed.x)*A + B}))
View(xrev2)
