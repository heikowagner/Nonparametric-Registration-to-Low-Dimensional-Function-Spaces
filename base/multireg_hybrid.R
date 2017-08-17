###Registration only over the functional space via optimize
#require(fda)
#library(trust)

require(orthogonalsplinebasis)
require(minqa)

#Xn - Unwarped Data
#K - Dimension of the warped space, if left out the dimension is estimated via algorithm 1
#w2 - Scaling paramter 
#d - abort criteria for algorithm 1, if not specified estiamted error variance
#T2 - Number of datapoint used for linear interpolation
#iter - Number of iterations performed by newoua algorithm
#B - Number of Basisfunction used to approximate W


getlam2 = function(scorew, L) {
    scorew2 = c(matrix(scorew, nrow = B) - rowMeans(matrix(scorew, nrow = B)))
    scorew2 = matrix(scorew2, nrow = B) / X2weights
    scorew2 = t(t(matrix(scorew2, nrow = B)) - rowMeans(t(matrix(scorew2, nrow = B))))
    hc = warph(scorew2)
    Xpre = sapply(1:N, function(i) Xf[[i]](hc[, i]))
    X = Xpre - Sp %*% Xpre;
    if (N < T) {
        s = svd(X, 0, 0)$d ^ 2
    } else {
        s = svd(t(X), 0, 0)$d ^ 2
    }
    out = ((1 - w2) * sum(s[ - (1:L)]) / sum(sg) + w2 / N * sum(scorew ^ 2))
    if (plotit == 1) {
        par(mfrow = c(1, 2))
        matplot(t2, X, typ = "l")
        matplot(t2, hc, typ = "l")
        print(out)
    }
    out
}

warph = function(scoresw) {
    ###later speedup due to : X2%*%t(ct)
    #h=sapply(1:N, function(i) cumsum(exp(X2%*%scoresw[((i-1)*B+1):(i*B)]))/sum(exp(X2%*%scoresw[((i-1)*B+1):(i*B)]))      )
    h = sapply(1:N, function(i) cumsum(exp(X2 %*% scoresw[((i - 1) * B + 1):(i * B)])))
    h = sapply(1:N, function(i)(h[, i] - min(h[, i])) / (max(h[, i] - min(h[, i]))))
}


multireg_hybrid = function(Xn, K, w2 = 0.01, d, T2 = 128, iter = 3000, B = 7, plotit = 1, c, basis) {

    B <<- B
    N <<- length(Xn);
    w2 <<- w2
    plotit <<- plotit
    Ba <<- basis
    Sp <<- Ba %*% solve(t(Ba) %*% Ba) %*% t(Ba)

    #grid of the original data
    #t=seq(from=0,to=1,length=T);

    #grid to approximate the data
    t2 <<- seq(from = 0, to = 1, length = T2);
    #Construct the orthonormal Basis functions for W
    interior = seq(from = 0, to = 1, length = (B - 2));
    knots = expand.knots(interior, order = 4);
    obase = OBasis(knots);
    X2 = evaluate(obase, t2);
    scaler = mean(diag(t(X2) %*% X2))
    X2 <<- X2 %*% diag(1 / (sqrt(diag(t(X2) %*% X2))))

    X2weights <<- colMeans(X2)
    #Estimate warping functions from c


    #Get linear interpolation of the data 
    #Xf=sapply(1:N, function(i) approxfun(Xn[[i]]$x,Xn[[i]]$y,rule = 2,method="linear"))
    Xf <<- sapply(1:N, function(i) approxfun((Xn[[i]]$x - min(Xn[[i]]$x)) / (max(Xn[[i]]$x) - min(Xn[[i]]$x)), Xn[[i]]$y, rule = 2, method = "linear"))

    #Estimate (w2-1) S(K,h(c)) + w2 V
    X = sapply(1:N, function(i) Xf[[i]](t2))
    sg <<- sum(svd(Sp %*% X, 0, 0)$d ^ 2)

    S_min = newuoa(c, getlam2, L = K, control = list(maxfun = iter, rhobeg = 0.001))

    scorew = S_min$par
    scorew2 = t(t(matrix(scorew, nrow = B)) - rowMeans(t(matrix(scorew, nrow = B))))
    scorew2 = c(matrix(scorew2, nrow = B) - rowMeans(matrix(scorew2, nrow = B)))
    c_star = matrix(scorew2, nrow = B) / X2weights
    h = warph(c_star)
    Xreg = sapply(1:N, function(i) Xf[[i]](h[, i]))

    list(h = h, X = Xreg, c = c_star, varW = var(S_min$par), Sstar = ((S_min$fval - w2 * var(S_min$par)) / (1 - w2)), c_pass = scorew)
}
