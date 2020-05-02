#--------------------------------------------------------------------------------
# File: RegressionHelpers.cs
# Author: Timothy O'Connor
# © Copyright University of Queensland, 2012-2014. All rights reserved.
# License: 
#--------------------------------------------------------------------------------

regress <- function(file, formula)
{
	x = read.table(file);
	xn = log(x + 1);

	print(dim(xn))

	r = lm(formula, data=xn);
	
	ra = summary(r);
	
	adjR <- sprintf("#0#%f#", ra[9])
	
	print(adjR);
}

#regress.data <- function(file, formula, columns, genes, geneSets)
#{
#	x = read.table(file);
#	xn = log(x + 1);
#
#	print(dim(xn))
#
#	r = lm(formula, data=xn);
#	
#	return r;
#}

permutation.test <- function(x1, y1, nshuffle)
{
	#DV <- c(x1, y1)
    #IV <- factor(rep(c("A", "B"), c(length(x1), length(y1))))
	#library(coin)   ;                 # for oneway_test(), pvalue()
	##return(pvalue(oneway_test(DV ~ IV, alternative="greater",
    #             distribution=approximate(B=999999))));
    
    meandif <- abs(mean(x1) - mean(y1));
    
    data <- c(x1, y1);
    
    s1 <- length(x1);
    s2 <- length(y1);
    
    
    #tst <- sum(boot(data, testfunc, nshuffle)$t);
    #print(sprintf("tstsum: %d; nshuffle: %d; pval %f == %e", tst, nshuffle, tst/nshuffle, tst/nshuffle))
    
    testfunc <- function(x, d) { s <- sample(s1+s2, s1+s2); sampmeandiff <- abs(mean(data[s[1:s1]]) - mean(data[s[(s1+1):(s1+s2)]])); if (sampmeandiff > meandif ) { return(1) } else { return(0) } }
    
    dummy <- vector(length=1);
    return(sum(boot(dummy, testfunc, nshuffle, parallel="multicore", ncpus=4)$t)/nshuffle)
}

regress.drm.boot <- function(fileMap, fileNN, fileNNR, fileDR, listSampleFiles, columnsTSS, columnsDRM, repCount, batchSize)
{
	library(boot)

	columns <- c(columnsTSS, columnsDRM);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	
	formulaTSS <- as.formula(paste("esc ~ ", paste(columnsTSS, collapse = "+" )))
    formulaDRM <- as.formula(paste("esc ~ ", paste(columnsDRM, collapse = "+" )))

	files <- c(fileMap, fileNN, fileNNR);
	
	if (length(listSampleFiles) > 0)
	{	files <- c(fileMap, fileNN, fileNNR, fileDR); }
	
	formulae <- list(formula, formulaTSS, formulaDRM);
	columnGroups <- list()
	columnGroups[[1]] = columns;
	columnGroups[[2]] = columnsTSS;
	columnGroups[[3]] = columnsDRM;
	
	load <- function(x) { return(read.table(x)) };
	#logd <- function(x) { return(log(x + 1)); };
	
	xd <- lapply(files, load);
	xs <- list();
	for (i in 1:length(xd)) { xs[[i]] <- log(xd[[i]] + 1) }
	
	#source(listScriptFile);
	#drmListData <- drmdata();
	
	size = dim(xs[[1]])[1]
	
	B <- list()
	
	for (i in 1:2)
	{
		B[[i]] <- list()
	
		for (j in 1:3)
		{
			form <- formulae[[j]];
			cols <- columnGroups[[j]]
			
			r <- lm(form, xs[[i]]);
			
			y <- as.matrix(xs[[i]]$esc);
			
			x <- as.matrix(xs[[i]][,cols])
			
			for (k in 1:length(r$coefficients[cols]))
			{
				if (is.na(r$coefficients[cols][k]))
				{
					r$coefficients[cols][k] = 0;
				}
			}
			
			res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[cols])
			tot <- y - mean(y)
			
			res2 <- res^2
			tot2 <- tot^2
			
			fstat <- function(x, d) { r2 <- 1 - sum(res2[d]) / sum(tot2[d]); adjr2 <- 1 - (1 - r2) * (length(y) - 1) / (length(y) - length(cols) - 1); return(adjr2);  }
		
			print(sprintf("i = %d; j = %d", i, j));
		
			if (i == 2 && j == 2)
			{
				# Don't repeat TSS only data
				B[[i]][[j]] <- B[[1]][[2]];
			}
			else
			{
				if (j != 3)
				{
					B[[i]][[j]] <- boot(res2, fstat, repCount)
				}
				else
				{
					dstat <- function(x, d) { return(0); }
					B[[i]][[j]] <- boot(res2, dstat, repCount)
				}
			}
		}
	}
	
	B[[3]] <- list()
	#rstat <- function(x, d) { return(summary(lm(formula, data=d[x[1:size],]))$adj.r.squared) };
	#rstatTSS <- function(x, d) { return(summary(lm(formulaTSS, data=d[x[1:size],]))$adj.r.squared) };
	#rstatDRM <- function(x, d) { return(summary(lm(formulaDRM, data=d[x[1:size],]))$adj.r.squared) };
	rstat <- function(x, d) { return(summary(lm(formula, data=xs[[3]][sample(dim(xs[[3]])[1],size),]))$adj.r.squared) };
	rstatTSS <- function(x, d) { return(summary(lm(formulaTSS, data=xs[[3]][sample(dim(xs[[3]])[1],size),]))$adj.r.squared) };
	rstatDRM <- function(x, d) { return(0); }; #return(summary(lm(formulaDRM, data=xs[[3]][sample(dim(xs[[3]])[1],size),]))$adj.r.squared) };
	
	dummy <- vector(length=1);
	B[[3]][[1]] <- boot(dummy, rstat, repCount * 10, parallel="multicore", ncpus=6);
	B[[3]][[2]] <- boot(dummy, rstatTSS, repCount * 10, parallel="multicore", ncpus=6);
	B[[3]][[3]] <- boot(dummy, rstatDRM, repCount * 10);
	
	
	B[[4]] <- list();
	
	#ns <- as.matrix(sequence(size));
	#fn <- function(x, y) { return ( sum ( drmListData[[x]][[y]][ sample( length(drmListData[[x]][[y]]), xd[[4]][x, "Count"] ) ] ) ) }
	#
	#sstat <- function(x, d) { 
	#    
	#    for (m in 1:length(columnsDRM)) { xs[[4]][,columnsDRM[m]] <- apply(ns, 1, fn, y=m) }
	#    
	#    return(summary(lm(formula, data=xs[[4]]))$adj.r.squared);
	#    };
	    
    #B[[4]][[1]] = boot(dummy, sstat, repCount);
	
	maxDataSets <- 3;
	if (length(listSampleFiles) > 0)
	{
		maxDataSets <- 4;
		B[[4]][[1]] <- list()
		B[[4]][[1]]$t <- matrix(NA, length(listSampleFiles) * batchSize, 1)
		setindex <- 1;
		for (i in 1:length(listSampleFiles))
		{
			x2 <- read.table(listSampleFiles[[i]]);
			
			xs2 <- log(x2 + 1);
			
			for (j in 1:batchSize)
			{
				B[[4]][[1]]$t[setindex,1] <- summary(lm(formula, xs2[(size * (j-1) + 1):(size * j),]))$adj.r.squared;
				setindex <- setindex + 1
			}
		}
		B[[4]][[1]]$t0 <- mean(B[[4]][[1]]$t)
		print(sprintf("Rsquareds|%d|%d|%f|", 4, 1, B[[4]][[1]]$t0));
		
	}
	else
	{
		print(sprintf("Rsquareds|%d|%d|%f|", 4, 1, 0));
	}
	
	for (i in 1:3) { for (k in 1:3) { print(sprintf("Rsquareds|%d|%d|%f|", i, k, B[[i]][[k]]$t0)); } }
	
	
	for (k in 1:3) { for (l in 1:3) { for (i in 1:maxDataSets) { for (j in 1:maxDataSets) {
		if ((i < 4 || k == 1) && (j < 4 || l == 1))
		{
			pval = 1.0;
			if (i != j || k != l)
			{
				pval = permutation.test(B[[i]][[k]]$t, B[[j]][[l]]$t, 100000);
			}
			print(sprintf("pvalue|%d|%d|%d|%d|%e|", i, k, j, l, pval));
		}}}
			
	}
	}
	
}

regress.drm <- function(file, formula, formulaTSS, formulaDRM)
{
	x = read.table(file);
	xn = log(x + 1);

	print(dim(xn))

	rstat <- function(x, d) { return(summary(lm(formula, data=xn[sample(dim(xn)[1],size),]))$adj.r.squared) };
	
	r = lm(formula, data=xn);
	rTSS = lm(formulaTSS, data=xn);
	rDRM = lm(formulaDRM, data=xn);

	ra = summary(r);
	raTSS = summary(rTSS);
	raDRM = summary(rDRM);

	adjR <- sprintf("#0#%f#", ra[9])
	adjRTSS <- sprintf("#1#%f#", raTSS[9])
	adjRDRM <- sprintf("#2#%f#", raDRM[9])

	#sink("~/tempROutput.txt", append=FALSE);
	
	print(adjR);
	print(adjRTSS);
	print(adjRDRM);
}

#regress.columns.r2sample <- function(file, columns, repCount, r2filename)
regress.columns.DUMMY <- function(file, columns, repCount, r2filename)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	r <- lm(formula, xn);
	
	print(sprintf("AdjR: %f", summary(r)$adj.r.squared));
			
	y <- as.matrix(xn$esc);
	x <- as.matrix(xn[,columns])
			
	for (k in 1:length(r$coefficients[columns]))
	{
		if (is.na(r$coefficients[columns][k]))
		{
			r$coefficients[columns][k] = 0;
		}
	}
			
	res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[columns])
	tot <- y - mean(y)
	
	res2 <- res^2
	tot2 <- tot^2
	
	fstat <- function(x, d) { s <- sample(length(res2), length(res2), replace=TRUE); r2 <- 1 - sum(res2[s]) / sum(tot2[s]); adjr2 <- 1 - (1 - r2) * (length(y) - 1) / (length(y) - length(columns) - 1); return(adjr2);  }

	dummy <- vector(length=1);
	r2s <- boot(dummy, fstat, repCount)$t;
	
	save(list="r2s", file=r2filename)
		
	for (i in 1:length(r2s))
	{
		print(sprintf("Rsquareds|%d|%f|", i, r2s[i]));
	}
}

regress.columns.loocv <- function(file, columns, r2filename)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	size <- dim(xn)[1];

	if (size > 1)
	{
		res2 <- vector(length=size);
		tot2 <- vector(length=size);

		sample <- 2:size

		for (i in 1:size)
		{
			xtest <- xn[i,];
			xtrain <- xn[sample,];

			r <- lm(formula, xtrain);
			 r$coefficients[1]
			print(sprintf("AdjR: %f, intercep: %f", summary(r)$adj.r.squared, r$coefficients[1]));
			
			y <- as.matrix(xtest$esc);
			x <- as.matrix(xtest[,columns])
			
			for (k in 1:length(r$coefficients[columns]))
			{
				if (is.na(r$coefficients[columns][k]))
				{
					r$coefficients[columns][k] = 0;
				}
			}
			
			res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[columns])
			tot <- y - mean(xtrain$esc)
	
			res2[i] <- res^2
			tot2[i] <- tot^2
		
			# move the current test index into the 
			#  training set for the next iteration
			if (i < size)
			{
				sample[i] <- i
			}
		}
	
		
		#r2test <- 1 - sum(res2) / sum(tot2);
		
		r2s <- res2 # sqrt(sum(res2)/size)

		r2s <- 1 - sum(res2)/sum(tot2)
			
		save(list="r2s", file=r2filename)
		
		for (i in 1:length(r2s))
		{
			print(sprintf("Rsquareds|%d|%f|", i, r2s[i]));
		}
	}
	else
	{
		print(sprintf("Rsquareds|0|-inf|"));
	}
}

regress.columns.cv <- function(file, columns, r2filename)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	size <- dim(xn)[1];

	folds <- 10
	foldsize <- floor(size/folds)

	shuffle <- sample(size)

	l <- list()
	for (i in 1:folds)
	{
		l[[i]] <- list()
	}

	for (i in 1:size)
	{
		bin <- floor((i-1)*folds/size)+1
		cat(bin)	

		l[[bin]] <- c(l[[bin]], shuffle[i])
	}

	if (folds > 1)
	{
		res2 <- vector(length=size);
		tot2 <- vector(length=size);

		sample <- 2:folds

		for (i in 1:folds)
		{
			testindices <- unlist(l[i])
			
			xtest <- xn[testindices,];
			xtrain <- xn[unlist(l[sample]),];

			r <- lm(formula, xtrain);
			 r$coefficients[1]
			print(sprintf("AdjR: %f, intercep: %f", summary(r)$adj.r.squared, r$coefficients[1]));
			
			y <- as.matrix(xtest$esc);
			x <- as.matrix(xtest[,columns])
			
			for (k in 1:length(r$coefficients[columns]))
			{
				if (is.na(r$coefficients[columns][k]))
				{
					r$coefficients[columns][k] = 0;
				}
			}
			
			res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[columns])
			tot <- y - mean(xtrain$esc)
	
			for (j in 1:length(testindices))
			{
				res2[testindices[j]] <- res[j]^2
				tot2[testindices[j]] <- tot[j]^2
			}
		
			# move the current test index into the 
			#  training set for the next iteration
			if (i < size)
			{
				sample[i] <- i
			}
		}
	
		
		#r2test <- 1 - sum(res2) / sum(tot2);
		
		r2s <- res2 # sqrt(sum(res2)/size)

		r2s <- 1 - sum(res2)/sum(tot2)
			
		save(list="r2s", file=r2filename)
		
		for (i in 1:length(r2s))
		{
			print(sprintf("Rsquareds|%d|%f|", i, r2s[i]));
		}
	}
	else
	{
		print(sprintf("Rsquareds|0|-inf|"));
	}
}

############################################
# Cross-validated Lasso regression R2 log10 + 1 transform #
############################################
regress.columns.lasso.cv <- function(file, columns, r2filename, plotfile)
{
	library(lars)
	#library(foreach)
	#library(doMC)
	#registerDoMC(6)

	x <- read.table(file);
	xn <- log10(x+1);
	
	xc <- xn[,columns];
	
	y <- xn['esc'];
	
	size <- dim(xn)[1];

	folds <- 10
	foldsize <- floor(size/folds)

	shuffle <- sample(size)

	l <- list()

	binsize <- size / folds;
	start <- 1;
	for (i in 1:folds)
	{
		end <- start + binsize
		startindex = (floor(start) + 1)
		endindex = min(length(shuffle), (floor(end)))
		l[[i]] <- shuffle[startindex:endindex];
		start <- end
	}

	if (folds > 1)
	{
		yhat <- vector(length=size);
		res2 <- vector(length=size);
		tot2 <- vector(length=size);

		r2sample <- vector(length=folds);

		#foreach (i=1:folds) %dopar%
		for (i in 1:folds)
		{
			print(sprintf("Fold %d", i));
			testindices <- unlist(l[i])
			
			xtest <- xc[testindices,];
			ytest <- y[testindices,];

			trainindices = unlist(l[-i]);
			xtrain <- xc[trainindices,];
			ytrain <- y[trainindices,];
			
			rc <- try(cv.lars(as.matrix(xtrain), as.matrix(ytrain), index=seq(0,1,length=(length(columns)+1)), K=10));

			sindex <- which.min(rc$cv)

			r <- lars(as.matrix(xtrain), as.matrix(ytrain), type='lasso');

			predictions <- predict.lars(r, as.matrix(xtest), sindex/(length(columns)+1), mode="fraction")$fit
			res <- ytest - predictions
			tot <- ytest - mean(ytest)
	
			r2sample[i] <- 1 - sum(res ^ 2)/sum(tot ^ 2)

			print(sprintf("Fold R2: %s", r2sample[i]));

			for (j in 1:length(testindices))
			{
				res2[testindices[j]] <- res[j]^2
				tot2[testindices[j]] <- tot[j]^2
				yhat[testindices[j]] <- predictions[j]
			}
		}
		
		#r2test <- 1 - sum(res2) / sum(tot2);
		#r2s <- res2 # sqrt(sum(res2)/size)
		#r2s <- 1 - sum(res2)/sum(tot2)
		#save(list="r2s", file=r2filename)

		pdf(file=plotfile, width=6, height=6);
		plot(yhat, y[,1]);
		dev.off();

		print(sprintf("Rsquareds|0|%f|", mean(r2sample)));
		print(sprintf("Error|0|%f|", sd(r2sample)/sqrt(folds)));
	}
	else
	{
		print(sprintf("Rsquareds|0|-inf|"));
		print(sprintf("Error|0|-inf|"));
	}
}

############################################
# Lasso regression R2 log10 + 1 transform #
############################################
regress.columns.lasso <- function(file, columns, r2filename, plotfile)
{
	library(boot);
	library(lars)
	x <- read.table(file);
	xn <- x
	xn$esc <- log10(xn$esc + 1)

	xc <- log10(xn[,columns] + 1);
	
	y <- xn[,'esc'];

	tries <- 0;
	isch <- TRUE;
	r <- 'a';
	while (tries < 50 && isch)
	{
		r <- try(cv.lars(as.matrix(xc), as.matrix(y), index=seq(0,1,length=(length(columns)+1)), K=10)); #length(as.matrix(y)));
		tries <- tries + 1;
		isch = is.character(r)
	}
	print(tries)

	# retry since error may be transient
	if (is.character(r))
	{
		print(sprintf("Rsquareds|%d|%f|", 1, 0));
		print(sprintf("Error|%d|%f|", 1, 0));
	}
	else
	{
		
		sindex <- which.min(r$cv)

		r2 <- 1 - r$cv[which.min(r$cv)]/var(y)
		error <- 1 - (r$cv[sindex] + r$cv.error[sindex])/var(y)

		#print("Hi");
		#rc <- try(lars(as.matrix(xc), as.matrix(y), type='lasso'));
		#
		#predictions <- predict.lars(rc, as.matrix(xc), sindex/(length(columns)+1), mode="fraction")$fit
		#
		#pdf(file=plotfile, width=6, height=6);
		#plot(predictions, y);
		#dev.off();

		print(sprintf("Rsquareds|%d|%f|", 1, r2));
		print(sprintf("Error|%d|%f|", 1, r2 - error));
	}
}

######################################
# Lasso regression R2 inv.logit(log) #
######################################
regress.columns.logit.lasso <- function(file, columns, r2filename)
{
	library(boot);
	library(lars)
	x <- read.table(file);
	xn <- inv.logit(log(as.matrix(x)))*2-1

	xc <- xn[,columns];

	for (i in 1:dim(xc)[2])
	{
		xc[which(is.nan(xc[,i])),i] <- 0
	}

	y <- xn[,'esc'];
	
	r <- try(cv.lars(as.matrix(xc), as.matrix(y), K=10)); #length(as.matrix(y)));
	
	if (is.character(r))
	{
		print(sprintf("Rsquareds|%d|%f|", 1, 0));
	}
	else
	{
		sindex <- which.min(r$cv)
		
		r2 <- 1 - r$cv[which.min(r$cv)]/var(y)

		print(sprintf("Rsquareds|%d|%f|", 1, r2));
	}
}

############################################
# Linear regression R2 log10 + 1 transform #
############################################
regress.columns.r2 <- function(file, columns, predictionFile)
{
	library(boot);
	x <- read.table(file);
	xn <- log10(x + 1);

	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	r <- lm(formula, as.data.frame(xn));

	print(sprintf("Rsquareds|%d|%f|", 1, summary(r)$adj.r.squared));

	p <- predict(r);

	write.csv(p, file=predictionFile);	
}

########################################################
# Linear regression R2 log10 + 1 transform predictions #
########################################################
regress.columns.predict <- function(file, columns, predictionFile)
{
	library(boot);
	x <- read.table(file);
	xn <- log10(x + 0.00001);

	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	r <- lm(formula, as.data.frame(xn));

	p <- predict(r);

	l <- length(columns);
	c <- vector(length=l);

	for (i in 1:l)
	{
		c[i] = cor(xn[,i], xn$esc);
	}

	write.csv(c, file='../temp/cors.csv');

	write.csv(p, file=predictionFile);
}

#######################################
# Linear regression R2 inv.logit(log) #
#######################################
regress.columns.logit.r2 <- function(file, columns, predictionFile)
{
	library(boot);
	x <- read.table(file);
	xn <- as.matrix(inv.logit(log(as.matrix(x)))*2-1);

	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	r <- lm(formula, as.data.frame(xn));

	print(sprintf("Rsquareds|%d|%f|", 1, summary(r)$adj.r.squared));	
}

regress.columns.shuffleExpression <- function(file, columns, r2filename)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	xn$esc <- as.matrix(xn$esc[sample(dim(xn)[1])]);

	r <- lm(formula, xn);
	
	print(sprintf("AdjR: %f", summary(r)$adj.r.squared));
			
	y <- as.matrix(xn$esc);
	x <- as.matrix(xn[,columns])
			
	for (k in 1:length(r$coefficients[columns]))
	{
		if (is.na(r$coefficients[columns][k]))
		{
			r$coefficients[columns][k] = 0;
		}
	}
			
	res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[columns])
	tot <- y - mean(y)
	
	res2 <- res^2
	tot2 <- tot^2
	
	r2 <- 1 - sum(res2) / sum(tot2);
	adjr2 <- 1 - (1 - r2) * (length(y) - 1) / (length(y) - length(columns) - 1); 


	if (adjr2 == Inf)
	{
		adjr2 = 10000;
	}
	print(sprintf("Rsquareds|%d|%f|", 1, summary(r)$adj.r.squared));	
}



regress.columns.residuals <- function(file, columns, r2filename)
#regress.columns.r2sample <- function(file, columns, r2filename)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	r <- lm(formula, xn);
	
	print(sprintf("AdjR: %f", summary(r)$adj.r.squared));
			
	y <- as.matrix(xn$esc);
	x <- as.matrix(xn[,columns])
			
	for (k in 1:length(r$coefficients[columns]))
	{
		if (is.na(r$coefficients[columns][k]))
		{
			r$coefficients[columns][k] = 0;
		}
	}
			
	res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[columns])
	tot <- y - mean(y)
	
	res2 <- res^2
	tot2 <- tot^2
	
	r2 <- 1 - sum(res2) / sum(tot2);
	adjr2 <- 1 - (1 - r2) * (length(y) - 1) / (length(y) - length(columns) - 1); 

	adjr2 <- res2
	r2s <- adjr2

	save(list="adjr2", file=r2filename)
	for (i in 1:length(r2s))
	{
		print(sprintf("Rsquareds|%d|%f|", i, r2s[i]));
	}
}

test.ranked.matrix <- function(filename)
{
	x <- read.table(filename)

	rows = dim(x)[1]
	cols = dim(x)[2]

	c <- matrix(0, nrow=rows,ncol=cols)
	for (i in 1:cols) { for (j in 1:rows) { if (x[j,1] > x[j,i]) { c[j,i] = c[j,i] + 1 } } }

	#for (i in 1:cols) { print(pbinom(max(sum(c[,i]),rows-sum(c[,i])),rows,0.5)) }

	#ct <- 0
	#for (i in 1:rows) { if (max(x[i,])==x[i,1]) { ct <- ct + 1 } }
	#print(pbinom(min(rows*cols - rows -sum(c), sum(c)), rows*cols-rows, 0.5))
	#print(1-pbinom(ct,rows,1/cols))

	for (t in 1:10) 
	{
		ct <- 0; 
		for (i in 1:rows) 
		{ 
			f <- which(x[i,] == max(x[i,]));
			if (f <= t) { ct <- ct + 1 }
		}
		p <- binom.test(ct,rows,t/cols,alternative="two.sided")$p.value 
		if (p == TRUE)
		{
			p <- 1;
		}
		print(p);
	}
}

regress.columns.rowsample <- function(file, columns, repCount, rowCount, r2filename)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	

	rstat <- function(x, d) { return(summary(lm(formula, xn[sample(dim(xn)[1],rowCount),]))$adj.r.squared) };
	
	dummy <- vector(length=1);
	if (.Platform$OS.type == "unix")
	{
		r2s <- boot(dummy, rstat, repCount)$t;
	}
	else
	{
		r2s <- boot(dummy, rstat, repCount, parallel="multicore", ncpus=6)$t;
	}

	save(list="r2s", file=r2filename)
		
	for (i in 1:length(r2s))
	{
		print(sprintf("Rsquareds|%d|%f|", i, r2s[i]));
	}
}

test.r2diff <- function(r2file1, r2file2)
{
    load(r2file1);
    r2null <- r2s;
    lead(r2file2);
    
    pval <- permutation.test(r2null, r2s, 100000);
    print(sprintf("pval|%f|", pval));
}

regress.promandmap <- function(file, columnsTSS, columnsDRM, repCount)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	columns <- c(columnsTSS, columnsDRM);
	
	formula <- as.formula(paste("esc ~ ", paste(columns, collapse = "+" )))	
	formulaTSS <- as.formula(paste("esc ~ ", paste(columnsTSS, collapse = "+" )))

	formulae <- list(formulaTSS, formula);
	columnGroups <- list(columnsTSS, columns)

	R2s <- list();
	
	i <- 1;
		form <- formulae[[i]];
		cols <- columnGroups[[i]]
			
		r <- lm(form, xn);
			
		y <- as.matrix(xn$esc);
		x <- as.matrix(xn[,cols])
			
		for (k in 1:length(r$coefficients[cols]))
		{
			if (is.na(r$coefficients[cols][k]))
			{
				r$coefficients[cols][k] = 0;
			}
		}
			
		res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[cols])
		tot <- y - mean(y)
		
		res2 <- res^2
		tot2 <- tot^2
		
		fstat <- function(x, d) { r2 <- 1 - sum(res2[d]) / sum(tot2[d]); adjr2 <- 1 - (1 - r2) * (length(y) - 1) / (length(y) - length(cols) - 1); return(adjr2);  }
	
		dummy <- vector(length=1);
		R2s[[i]] <- boot(dummy, fstat, repCount);
		
		print(sprintf("Rsquareds|%d|%f|", i, R2s[[i]]$t0));
	
	sstat <- function(x, d) 
	{
		l <- dim(xn)[1];
		data <- c(xn[,columnsTSS], xn[sample(l),columnsDRM], xn$esc)
		r2 <- summary(lm(formula, data))$adj.r.squared;
		rm(data);	
		gc();
		return (r2);
	}
	
	data <- xn;
	rvec <- vector(length=repCount);
	
	for (i in 1:repCount)
	{
		l <- dim(xn)[1];
		for (j in 1:length(columnsDRM))
		{
			data[,columnsDRM[j]] <- xn[sample(l),columnsDRM[j]]
		}
		#c(xn[,columnsTSS], xn[sample(l),columnsDRM], xn$esc)
		r2 <- summary(lm(formula, data))$adj.r.squared;
		#rm(data);	
		gc();
	
		rvec[i] <- r2;
	}
	#R2s[[2]] <- boot(dummy, sstat, repCount);
	
	pval <- permutation.test(R2s[[1]]$t, R2s[[2]]$t, 10000);
	print(sprintf("pvalue|%e|", pval));
}

regress.promandmap2 <- function(file, columnsTSS, columnsDRM, columnsDRM2, repCount)
{
	library(boot);
	x <- read.table(file);
	xn <- log(x + 1);
	
	columns1 <- c(columnsTSS, columnsDRM);
	columns2 <- c(columnsTSS, columnsDRM2);
	
	formula1 <- as.formula(paste("esc ~ ", paste(columns1, collapse = "+" )))	
	formula2 <- as.formula(paste("esc ~ ", paste(columns2, collapse = "+" )))	
	formulaTSS <- as.formula(paste("esc ~ ", paste(columnsTSS, collapse = "+" )))

	formulae <- list(formulaTSS, formula1, formula2);
	columnGroups <- list(columnsTSS, columns1, columns2)

	R2s <- list();
	
	for (i in 1:3)
	{
		form <- formulae[[i]];
		cols <- columnGroups[[i]]
			
		r <- lm(form, xn);
			
		y <- as.matrix(xn$esc);
		x <- as.matrix(xn[,cols])
			
		for (k in 1:length(r$coefficients[cols]))
		{
			if (is.na(r$coefficients[cols][k]))
			{
				r$coefficients[cols][k] = 0;
			}
		}
			
		res <- y - r$coefficients[1] - x %*% as.matrix(r$coefficients[cols])
		tot <- y - mean(y)
		
		res2 <- res^2
		tot2 <- tot^2
		
		fstat <- function(x, d) { r2 <- 1 - sum(res2[d]) / sum(tot2[d]); adjr2 <- 1 - (1 - r2) * (length(y) - 1) / (length(y) - length(cols) - 1); return(adjr2);  }
	
		R2s[[i]] <- boot(res2, fstat, repCount);
		
		print(sprintf("Rsquareds|%d|%f|", i, R2s[[i]]$t0));
	}
	pval <- permutation.test(R2s[[1]]$t, R2s[[2]]$t, 100000);
	print(sprintf("pvalue|%e|", pval));
	pval <- permutation.test(R2s[[1]]$t, R2s[[3]]$t, 100000);
	print(sprintf("pvalue|%e|", pval));
	pval <- permutation.test(R2s[[2]]$t, R2s[[3]]$t, 100000);
	print(sprintf("pvalue|%e|", pval));
}

regress.pcr.drm <- function(file, columns, columnsTSS, columnsDRM)
{
	library(plsdof)
	
	x = read.table(file);
	xn = log(x + 1);

	print(dim(xn))

	print(dim(xn))

	n <- dim(xn)[1];

    y <- as.vector(xn[,"esc"])

	xm <- as.matrix(xn[columns])
	xmTSS <- as.matrix(xn[columnsTSS])
	xmDRM <- as.matrix(xn[columnsDRM])

	r <- pcr(xm, y)
	rTSS <- pcr(xmTSS, y)
	rDRM <- pcr(xmDRM, y)

	#ra = summary(r);
	#raTSS = summary(rTSS);
	#raDRM = summary(rDRM);
	
	m <- length(columns);
	mTSS <- length(columnsTSS);
	mDRM <- length(columnsDRM);

	res <- matrix(, length(y), m + 1)
	resTSS <- matrix(, length(y), mTSS + 1)
	resDRM <- matrix(, length(y), mDRM + 1)

	for (j in 1:(m + 1)) {
		res[, j] <- y - r$intercept[j] - xm %*% r$coefficients[,j]
	}
	res2 <- apply(res^2, 2, sum)

	for (j in 1:(mTSS + 1)) {
		resTSS[, j] <- y - rTSS$intercept[j] - xmTSS %*% rTSS$coefficients[,j]
	}
	resTSS2 <- apply(resTSS^2, 2, sum)

	for (j in 1:(mDRM + 1)) {
		resDRM[, j] <- y - rDRM$intercept[j] - xmDRM %*% rDRM$coefficients[,j]
	}
	resDRM2 <- apply(resDRM^2, 2, sum)
	
	SStot <- sum((y - mean(y))^2)
	
	adjR <- sprintf("#0#%f#", 1 - res2[min(length(columns), 10)]/SStot)
	adjRTSS <- sprintf("#1#%f#", 1 - resTSS2[min(length(columnsTSS), 10)]/SStot)
	adjRDRM <- sprintf("#2#%f#", 1 - resDRM2[min(length(columnsDRM), 10)]/SStot)

	#sink("~/tempROutput.txt", append=FALSE);
	
	print(adjR);
	print(adjRTSS);
	print(adjRDRM);
}

regress.pcr.cv.drm <- function(file, columns, columnsTSS, columnsDRM)
{
	library(plsdof)

	x <- read.table(file);
	xn <- log(x + 1);

	print(dim(xn))

	n <- dim(xn)[1];

    y <- as.vector(xn[,"esc"])

	xm <- as.matrix(xn[columns])
	xmTSS <- as.matrix(xn[columnsTSS])
	xmDRM <- as.matrix(xn[columnsDRM])

	nFold <- 10

	r <- uq.pcr.cv(xm, y, k=nFold)
	rTSS <- uq.pcr.cv(xmTSS, y, k=nFold)
	rDRM <- uq.pcr.cv(xmDRM, y, k=nFold)

	SStot <- sum((y - mean(y))^2)

        
	adjR <- sprintf("#0#%f#%d#", 1 - r$cv.error[r$m.opt] * nFold / SStot, r$m.opt)
	if (r$m.opt == 0)
	{
		adjR <- sprintf("#0#0.0#0#");
	}
	adjRTSS <- sprintf("#1#%f#%d#", 1 - rTSS$cv.error[rTSS$m.opt] * nFold / SStot, rTSS$m.opt)
	if (rTSS$m.opt == 0)
	{
		adjRTSS <- sprintf("#1#0.0#0#");
	}
	adjRDRM <- sprintf("#2#%f#%d#", 1 - rDRM$cv.error[rDRM$m.opt] * nFold / SStot, rDRM$m.opt)
	if (rDRM$m.opt == 0)
	{
		adjRDRM <- sprintf("#2#0.0#0#");
	}

	#sink("~/tempROutput.txt", append=FALSE);
	
	print(adjR);
	print(adjRTSS);
	print(adjRDRM);
}

uq.pcr.cv <- function (X, y, k = 10, m = min(ncol(X), nrow(X) - 1), groups = NULL, 
    scale = TRUE, eps = 1e-06, plot.it = FALSE, compute.jackknife = TRUE) 
{
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(groups) == FALSE) {
        f = as.factor(groups)
        k = length(levels(f))
        my.names = levels(f)
    }
    if (is.null(groups) == TRUE) {
        f <- rep(1:k, length = n)
        my.names <- 1:k
    }
    all.folds <- split(sample(1:n), f)
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1)
    cv.error.matrix = matrix(0, k, m + 1)
    rownames(cv.error.matrix) = my.names
    colnames(cv.error.matrix) = 0:m
    coefficients.jackknife = NULL
    if (compute.jackknife == TRUE) {
        coefficients.jackknife <- array(dim = c(p, m + 1, k))
    }
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pcr.object <- pcr(Xtrain, ytrain, scale, m, eps)
        res <- matrix(, length(ytest), m + 1)
        if (compute.jackknife == TRUE) {
            coefficients.jackknife[, , i] = pcr.object$coefficients
        }
        for (j in 1:(m + 1)) {
            res[, j] <- ytest - pcr.object$intercept[j] - Xtest %*% 
                pcr.object$coefficients[, j]
        }
        cv.error.matrix[i, ] <- apply(res^2, 2, sum)
    }
    cv.error <- apply(cv.error.matrix, 2, mean)
    m.opt <- which.min(cv.error) - 1
    if (plot.it == TRUE) {
        plot(0:m, cv.error, type = "l")
    }
    mFull <- m.opt
    if (m.opt < 1)
    {
        mFull <- 1
    }
    pcr.object <- pcr(X, y, scale, mFull, eps = eps)
    coefficients <- pcr.object$coefficients[, m.opt + 1]
    intercept <- pcr.object$intercept[m.opt + 1]
    return(list(intercept = intercept, coefficients = coefficients, 
        m.opt = m.opt, cv.error.matrix = cv.error.matrix, cv.error = cv.error, 
        coefficients.jackknife = coefficients.jackknife))

	pcr.object <- pcr(X, y)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

