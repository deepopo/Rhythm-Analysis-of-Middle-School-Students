
cosinor <- function(t, y, w, alpha){
	n <- length(t);
	x <- cos(w * t);
	z <- sin(w * t);
	NE <- t(array(c(n, sum(x), sum(z), 
      sum(x), sum(x^2), sum(x*z), 
      sum(z), sum(x*z), sum(z^2)), dim=c(3, 3)))
	NEy = c(sum(y), sum(x*y), sum(z*y))
	###对以上矩阵进行化简，化简成行阶梯矩阵###
	RNE <- solve(NE, NEy)
	M <- RNE[1]; beta <- RNE[2]; gamma <- RNE[3];
	Amp <- sqrt(beta^2 + gamma^2);
	theta <- atan(abs(gamma / beta));
	a <- sign(beta);
	b <- sign(gamma);
	if ((a == 1 || a == 0) && b == 1)
	{
		phi <- -theta;
	}
	else if (a == -1 && (b == 1 || b == 0))
	{
		phi <- -pi + theta;
	}
	else if ((a == -1 || a == 0) && b == -1)
	{
		phi <- -pi - theta;
	}
	else if (a == 1 && (b == -1 || b == 0))
	{
		phi <- -2*pi + theta;
	}
	f <- M + Amp * cos(w * t + phi);
	R2 <- sum((f - mean(y)) ^ 2)/sum((y - mean(y)) ^ 2);
	Acr <- -phi * 24 / w / 2880;
	return(c(M, Amp, Acr, R2));
}
space <- function(x){
	return (paste(substr(x, 1, 13), substr(x, 14, 100)));
}
nroot <- function(x, n){
    return (abs(x) ^ n * sign(x));
}
ht <- function(ct, m, gama){
	return (nroot(ct, gama) / (nroot(m, gama) + nroot(ct, gama)));
}
lt <- function(ct, alpha, beta){
	return (exp(beta * (ct - alpha)) / (1 + exp(beta * (ct - alpha))));
}
xt <- function(ct, alpha, beta){
	return (atan(beta * (ct - alpha)) /pi + 0.5);
}
clip <- function(x){
	if (x < -1)
	{
		return (-1);
	} else if (x > 1)
	{
		return (1);
	} else{
		return (x);
	}
}
###加载颜色库###
library(RColorBrewer);
###需要多少种颜色###
cols<-brewer.pal(n=8,name="Set2");
###设置路径###
path <- 'XXXXXX/project/';
library(ggplot2);

w <- 2*pi*30/86400;
alpha <- 0.05;
###读取数据库###
files <- read.csv(paste(path, 'data/data.csv', sep=''), fileEncoding='UTF-8');
file_list <- list.files(path = paste(path, 'data', sep = ''), full.names=FALSE, recursive=TRUE)
file_name <- c(0);
Mesor <- c(0);
Amplitude <- c(0);
Acrophase <- c(0);
R2 <- c(0);
R2_h <- c(0);
R2_l <- c(0);
R2_x <- c(0);
amp_h <- c(0);
Min <- c(0);
xmin <- c(0);
m <- c(0);
gamma <- c(0);
t05u_h <- c(0);
t05l_h <- c(0);
width_ratio_h <- c(0)
amp_l <- c(0);
alpha_l <- c(0);
beta_l <- c(0);
t05u_l <- c(0);
t05l_l <- c(0);
width_ratio_l <- c(0);
amp_x <- c(0);
alpha_x <- c(0);
beta_x <- c(0);
t05u_x <- c(0);
t05l_x <- c(0);
width_ratio_x <- c(0);
len <- length(t(files))/length(files);
FILE <- 'MOS2D26170041 (2017-12-05)30sec.csv';
i <- 1; 

file_name[i] <- files[i, 1];
file <- space(files[i, 1]);
if (FILE %in% file_list){
	file <- read.csv(paste(path, 'R/', FILE, sep = ''));
	file <- file[!is.na(file[ , 2]), 2];
	file <- log(file + 1);
	d_num <- length(file) / 2880;
	###一天的时间点###
	t <- seq(2880 * d_num);
	x <- cosinor(t, file, w, alpha);
	Mesor[i] <- x[1];
	Amplitude[i] <- x[2];
	Acrophase[i] <- x[3];
	R2[i] <- x[4];
	cosine <- x[1] + x[2] * cos(w * t - x[3] * w * 2880 /24);
	ct <- cos(w * t - x[3] * w * 2880 /24);
	Min[i] <- x[1] - x[2];
	xmin[i] <- x[1] - 1.5 * x[2];

	model_ht <- try(nls(cosine ~ Min[[i]] + 2.8 * ht(cosine, m, gama), start = list(gama = 1.4, m = 0.5), control = nls.control(tol = 1e-03), algorithm = 'port'), silent = TRUE);
	if ('try-error' %in% class(model_ht)){
		m[i] <- NA; gamma[i] <- NA;
		t05u_h[i] <- NA;
		t05l_h[i] <- NA;
		width_ratio_h[i] <- NA;
		hill <- Min[[i]] + 2.8 * x[2] * ht(cosine, m[[i]], gamma[[i]]);
		R2_h[i] <- NA;
	} else{
		confi_ht <- summary(model_ht)$coefficients[ , 1];
		gamma[i] <- confi_ht[1]; m[i] <- confi_ht[2];
		t05u_h[i] <- x[3] + acos(clip(m[[i]] - x[2] + 0.5)) * (2 * pi / 24);
		t05l_h[i] <- x[3] - acos(clip(m[[i]] - x[2] + 0.5)) * (2 * pi / 24);
		width_ratio_h[i] <- (t05u_h[i] - t05l_h[i]) / 24;
		hill <- Min[[i]] + 2.8 * x[2] * ht(cosine, m[[i]], gamma[[i]]);
		R2_h[i] <- sum((hill - mean(file)) ^ 2)/sum((file - mean(file)) ^ 2);
	}
	model_lt <- try(nls(cosine ~ Min[[i]] + 2 * lt(ct, alpha, beta), start = list(alpha = 0, beta = 2), control = nls.control(tol = 1e-03), algorithm = 'port'), silent = TRUE);
	if ('try-error' %in% class(model_lt)){
		alpha_l[i] <- NA; beta_l[i] <- NA;
		t05u_l[i] <- NA;
		t05l_l[i] <- NA;
		width_ratio_l[i] <- NA;
		R2_l[i] <- NA;
	} else{
		confi_lt <- summary(model_lt)$coefficients[ , 1];
		alpha_l[i] <- confi_lt[1]; beta_l[i] <- confi_lt[2];
		t05u_l[i] <- x[3] + acos(clip(alpha_l[i])) * (2 * pi / 24);
		t05l_l[i] <- x[3] - acos(clip(alpha_l[i])) * (2 * pi / 24);
		width_ratio_l[i] <- (t05u_l[i] - t05l_l[i]) / 24;
		logistic <- Min[[i]] + 2 * x[2] * lt(ct, alpha_l[[i]], beta_l[[i]]);
		R2_l[i] <- sum((logistic - mean(file)) ^ 2)/sum((file - mean(file)) ^ 2);
	}
	model_xt <- try(nls(cosine ~ xmin[[i]] + 2 * xt(ct, alpha, beta), start = list(alpha = 0, beta = 2), control = nls.control(tol = 1e-03), algorithm = 'port'), silent = TRUE);
	if ('try-error' %in% class(model_xt)){
		alpha_x[i] <- NA; beta_x[i] <- NA;
		t05u_x[i] <- NA;
		t05l_x[i] <- NA;
		width_ratio_x[i] <- NA;
		R2_x[i] <- NA;
	} else{
		confi_xt <- summary(model_xt)$coefficients[ , 1];
		alpha_x[i] <- confi_xt[1]; beta_x[i] <- confi_xt[2];
		t05u_x[i] <- x[3] + acos(clip(alpha_x[i])) * (2 * pi / 24);
		t05l_x[i] <- x[3] - acos(clip(alpha_x[i])) * (2 * pi / 24);
		width_ratio_x[i] <- (t05u_x[i] - t05l_x[i]) / 24;
		arctangent <- xmin[[i]] + 2 * x[2] * xt(ct, alpha_x[[i]], beta_x[[i]]);
		R2_x[i] <- sum((arctangent - mean(file)) ^ 2)/sum((file - mean(file)) ^ 2);
	}
	#print(confi_xt);
	plot1f <- data.frame(x = t, y = file);
	plot2f <- data.frame(x = t, y = cosine);
	plot3f <- data.frame(x = t, y = hill);
	plot4f <- data.frame(x = t, y = logistic);
	plot5f <- data.frame(x = t, y = arctangent);

	plot <- ggplot()+
		geom_point(data = plot1f, aes(x = x, y = y, color = 'origin'), alpha = .2)+
		geom_line(data = plot2f, aes(x = x, y = y, color = 'cosine'), size = 1)+
		geom_line(data = plot3f, aes(x = x, y = y, color = 'hill'), size = 1)+
		geom_line(data = plot4f, aes(x = x, y = y, color = 'logistic'), size = 1)+
		geom_line(data = plot5f, aes(x = x, y = y, color = 'arctangent'), size = 1)+

	labs(title = 'Activity', y = 'activity', x = 'time');

		
}


df <- data.frame(file_name = file_name, Mesor = Mesor, Amplitude = Amplitude, Acrophase = Acrophase, R2 = R2, 
min_h = Min, amp_h = 2.8 * Amplitude, phi_h = Acrophase, m = m, gamma = gamma, t05u_h = t05u_h, t05l_h = t05l_h, width_ratio_h = width_ratio_h, R2_h = R2_h,
min_l = Min, amp_l = 2* Amplitude, phi_l = Acrophase, alpha_l = alpha_l, beta_l = beta_l, t05u_l = t05u_l, t05l_l = t05l_l, width_ratio_l = width_ratio_l, R2_l = R2_l,
min_x = xmin, amp_x = 2* Amplitude, phi_x = Acrophase, alpha_x = alpha_x, beta_x = beta_x, t05u_x = t05u_x, t05l_x = t05l_x, width_ratio_x = width_ratio_x, R2_x = R2_x)
print(df);
#write.table(df, paste(path, 'result/Parameter', '.csv', sep=''), row.names =  FALSE, col.names = c('file_name', 'Mesor', 'Amplitude', 'Acrophase', 'R2', 'min_h', 'amp_h', 'phi_h', 'm', 'gamma', 't05u_h', 't05l_h', 'width_ratio_h', 'R2_h', 'min_l', 'amp_l', 'phi_l', 'alpha_l', 'beta_l', 't05u_l', 't05l_l', 'width_ratio_l', 'R2_l', 'min_x', 'amp_x', 'phi_x', 'alpha_x', 'beta_x', 't05u_x', 't05l_x', 'width_ratio_x', 'R2_x'), sep = ',');

