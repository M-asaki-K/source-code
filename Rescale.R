#Re-scale of optimized x
A <- as.matrix(t(apply(x[ , ], 2, sd)))
B <- as.matrix(t(apply(x[ , ], 2, mean)))
#preprocessed.x�̕����͋K�i�����ꂽ�C�ӂ̃f�[�^�Z�b�g�ɕύX�\�i�Ⴆ�΋t��͂ŏo���œK�����Ƃ��j
xrev <- t(apply(preprocessed.x, 1, function(preprocessed.x){(preprocessed.x)*A + B}))
View(xrev2)