
# Allows the user to visually inspect the effect of injecting sigma noise to a normal dataset and
# compare its correlations to the real correlations on a tumor dataset. 

difconet.catf <- function(dObj, ...) {
	if (dObj$verbose) {
		cat(...)
		if (!is.null(flush.console)) flush.console()
	}
}


########################################################
#Function collection to create an artificial dataset   #
########################################################

# add noise to "normal" data (ndata) to compare with tumor data tdata
difconet.noise.inspection <- function(ndata, tdata, sigma=c(0.5, 0.75, 1.25), maxgenes=5000, corfunc=function(a,b) cor(a,b,method="spearman"))
{
	#corfunc=function(a,b) cor(a,b,method="spearman",use="pairwise.complete.obs")

	#ndata <- quantile.normalization(ndata)
	#tdata <- quantile.normalization(tdata)

	if (nrow(ndata) != nrow(tdata)) stop("ndata and tdata should contain the same number of rows.")

	if (nrow(ndata) > maxgenes) {
		wgenes <- sample(nrow(ndata), maxgenes)
		ndata <- ndata[wgenes, ]
		tdata <- tdata[wgenes, ]
	}

	N <- t(ndata)
	T <- t(tdata)
	rm(ndata)
	rm(tdata)

	#Calculate correlations
	N.cor <- corfunc(N, N)
	ndens <- density(N.cor[lower.tri(N.cor)], na.rm=TRUE)
	rm(N.cor)

	T.cor <- corfunc(T, T)
	tdens <- density(T.cor[lower.tri(T.cor)], na.rm=TRUE)
	rm(T.cor)

	plot(ndens, xlab="", main="Correlation Distributions", xlim=c(-1,1), ylim=c(0,max(ndens$y, tdens$y)), lwd=3)
	abline(v=0, col=8, lty=3)
	lines(tdens, col=8, lwd=2)

	for (i in 1:length(sigma)) {
		N.inj <- N +  matrix(rnorm(nrow(N)*ncol(N), mean=0, sd=sigma[i]), ncol=ncol(N))
		N.inj.cor <- corfunc(N.inj, N.inj)
		ninjdens <- density(N.inj.cor[lower.tri(N.inj.cor)], na.rm=TRUE)
		rm(N.inj.cor)
		lines(ninjdens, col=i+1, lwd=1, lty=2)
	}
	legend("topleft", c("Normal","Tumor",paste0("N+s=",sigma)), col=c(1,8,1:length(sigma)+1), lty=c(1,1,rep(2,length(sigma))), lwd=c(3,2,rep(1,length(sigma))), bg=0)
}


difconet.build.controlled.dataset <- 	function(
		data,
		noise.genes = round(nrow(data)*0.1),
		noise.sigma = c(0.0, 0.1, 0.2), # this is cummulative
		nonoise.sigma = c(0.0, 0.01, 0.01), # this is cummulative
		netcov = matrix(c(
			0.90, 0.90, 0.75, 0.75, 0.60, 0.60, 0.45, 0.45, 0.30, 0.30, 0.15, 0.15, 0.30, 0.30, 0.45, 0.45, 0.60, 0.60, 0.75, 0.75,
			0.95, 0.95, 0.80, 0.80, 0.65, 0.65, 0.50, 0.50, 0.35, 0.35, 0.10, 0.10, 0.25, 0.25, 0.40, 0.40, 0.55, 0.55, 0.70, 0.70,
			1.00, 1.00, 0.85, 0.85, 0.70, 0.70, 0.55, 0.55, 0.40, 0.40, 0.05, 0.05, 0.20, 0.20, 0.35, 0.35, 0.50, 0.50, 0.65, 0.65
			), ncol=3),
		genes.nets = 10,
		corfunc=function(a,b) cor(a,b,method="spearman"),
		verbose = TRUE
	)
{	
	sigma <- noise.sigma
	if(is.null(sigma)) {
		stop("A sigma value is needed. Run difconet.noise.inspection  to adjust noise level manually.")
	}
	if(ncol(netcov) != length(sigma)) {
		stop("Number of netcov column in networks has to be equal to the number of stages.")
	}
	if(length(genes.nets) > 1 &&  length(genes.nets) != nrow(netcov)) {
		stop("Number of genes in each net has to be equal to the number of networks (rows in netcov).")
	}
	if (length(genes.nets) == 1) genes.nets <- rep(genes.nets, nrow(netcov))
	N <- data

	if (is.null(rownames(N))) rownames(N) <- 1:nrow(N)
	
	rownames(N) <- paste0("Control_",rownames(N))

	rgenes <- c()
	if (noise.genes > 0) {
		rgenes <- sort(sample(nrow(N), noise.genes))
		if (verbose) cat(paste0("#Noisy Genes: " , length(rgenes), " out of ", nrow(N), "\n"))
		rn <- rownames(N)
		rn[rgenes] <- paste0("Noised_",rn[rgenes])
		rownames(N) <- rn
	}

	netidx <- nrow(N)
	if (sum(genes.nets) > 0) {
		mx <- matrix(0, ncol=ncol(N), nrow=sum(genes.nets))
		rownames(mx) <- paste0("Network_",1:nrow(mx))
		#netg <- 1:nrow(mx) + nrow(N)
		xN <- rbind(N, mx)
	} else {
		xN <- N+0
	}
	#require(mvtnorm)

	stages <-list()
	wno <- !(1:nrow(xN) %in% rgenes)

	for (i in 1:length(sigma)) {
		# Noise Injection
		if (noise.genes > 0 && sigma[i] > 0) {
			if (verbose) cat("Adding Noise To",length(rgenes),"genes for Stage", i," m=0, sd=",sigma[i],"\n")
			xN[rgenes, ] <- xN[rgenes, ] + matrix(rnorm(length(rgenes)*ncol(xN), mean=0, sd=sigma[i]), ncol=ncol(xN))
		}
		if (nonoise.sigma[i] > 0) {
			if (verbose) cat("Adding Noise To",sum(wno),"genes for Stage", i," m=0, sd=",nonoise.sigma[i],"\n")
			xN[wno, ] <- xN[wno, ] + matrix(rnorm(sum(wno)*ncol(xN), mean=0, sd=nonoise.sigma[i]), ncol=ncol(xN))
		}
		# Network Injection
		xI <- netidx
		for (j in 1:nrow(netcov)) {
			if (genes.nets[j] > 0) {
				sigma <- matrix(netcov[j,i], nrow=genes.nets[j], ncol=genes.nets[j])
				diag(sigma) <- 1
				if (verbose) cat("Generating Network Data", j, "Stage", i, ", Full Connected Networks of Correlation:", netcov[j,i], "\n")
				xM <- rmvnorm(ncol(xN), mean=rep(0,genes.nets[j]), sigma=sigma)
				#print(xM)
				xN[xI+1:genes.nets[j], ] <- t(xM)
				xI <- xI+genes.nets[j]
				if (verbose) {
					if (i==1 && j==1) {
						plot(density(corfunc(xN[1:xI,], xN[1:xI,]), na.rm=TRUE), col=1, xlim=c(-1,1))
					} else {
						lines(density(corfunc(xN[1:xI,], xN[1:xI,]), na.rm=TRUE), col=(i-1)*nrow(netcov)+j+1)
					}
					lines(density(corfunc(xM, xM), na.rm=TRUE), col=(i-1)*nrow(netcov)+j+1, lty=3)
				}
			}
		}
		# Store
		stages[[i]] <- xN
	}
	return (stages)
}



########################################################
#Function collection to analyse correlation differences#
########################################################

difconet.run <- function(	
						data,							#Data matrix : rownames with gene IDs colnames samples ID)
						predictor,						#Factor specifying the sample classes
						metric = c(1,2,3,4,5,6),		#Metrics to be calculated, if "8" is included, then a function is needed
						cutoff = 0.3, 					#Threshold for M1 and M3
						blocs = 5000,					#Dimension of the blocs when calculating the correlation matrices
						num_perms = 10,                 #Number of permutations to calculate a p-value
						comparisons = "all",           	#Must be a list with stage pairs to compare. e.g.: list(c(1,2), c(2,3)), or "all" for all possible 
														#combinations , or "seq" for sequential combinations. i.e.: 1-2,2-3,....,(n-1)-n
						perm_mode = "columns",          #How to perform the permutations: tags, tags by gene, 
						use_all_perm = TRUE,			#How the p value will be estimated, using all permutated data or only the ROW permutations (it requires more permutations)
						save_perm = FALSE,				#Save all permutated results?
						speedup = 0,					#Whether the speed-up mechanic should be used or not, 0 not speed up, 1:6 use this metric for speedup
						verbose = TRUE, 				#Whether info of the process in needed.
						metricfunc = NULL, 				#Function to estimate the distance metric needs a value "8"
						corfunc = function(a,b) cor(a,b,method="spearman")
						) 
{                             
	

	labels_vector <- as.character(predictor) #as.numeric(as.factor(predictor))	
	classes <- (unique(labels_vector))
	if (is.numeric(classes)) { classes <- sort(classes); }

	stagesL <- split.datset(classes, data, labels_vector, comparisons)
	stages <- stagesL[[1]]
	combinations <- stagesL[[2]]
	
	#Free up memory
	rm(data)
	gc()

	dObj <- list(
		combdens = list(),
		stages.data = stages,
		stage = stagesL[[3]],
		labels = classes,
		comparisons = comparisons,
		combinations = combinations,
		num_perms = num_perms,
		perm_mode = perm_mode,
		use_all_perm = use_all_perm,
		speedup = speedup,
		combstats = list(),
		verbose = verbose,
		corfunc = corfunc,
		metricfunc = metricfunc
					)
	
	#Start the process for each of the requested combinations

	for (i in 1:ncol(combinations)) {
	
		a <- stages[[ combinations[1,i] ]] #data.matrix(stages[[combinations[1,i]]])
		b <- stages[[ combinations[2,i] ]] #data.matrix(stages[[combinations[2,i]]])
		
		a_cols <- ncol(a)
		b_cols <- ncol(b)
		both_cols <- a_cols + b_cols
		
		#index  is a string representing the stages currently being analyzed (e.g.: 1-2)
		index <- paste(combinations[1,i], "_" ,combinations[2,i], sep = "")

		difconet.catf(dObj, "#################### COMBINATION ##################\n")
		difconet.catf(dObj, "Classes/Stages:",index,", samples(",combinations[1,i],"):",ncol(a),", samples(",combinations[2,i],"):",ncol(b),"\n")
	
		#M1, M2, M3, M3.1, M3.2..., M4
		metrics_and_cutoffs <- get.metric.cutoffs(dObj, metric, cutoff)
		#createDirs(metrics_and_cutoffs)			

		#Perform this only if a metric different from M2 is requested
		if(length(metric) > 0)
		{
			difconet.catf(dObj, "\nCalculating requested metrics...\n")
			#Gets the metric distance for the given correlation matrices a and b
			distance <- difconet.get.cor.distance(dObj, a, b, metric , cutoff, blocs)	
			#loss.gain <- difconet.get.cor.distance(a, b, "loss" , cutoff, blocs, full = TRUE)
			
			if(speedup)
			{
				difconet.catf(dObj, "\nDetermining Quantiles...\n")			
				#Split the distances according to their quantile.
				dist_quantiles <- apply(distance, 2, function(x) quantile(x,0:10/10, na.rm=TRUE))
				
				rep_metric <- which(metric == speedup)
				difconet.catf(dObj, "\nQuantile Sampling...\n")
				gene_samples <- get.representative.genes(distance[,rep_metric], dist_quantiles[,rep_metric]) 
				
				
				difconet.catf(dObj, "\nPermutations on sampled genes...\n")
				
				#Gets the permutations of the sampled genes for all the requested metrics
				samples_perm <- get.permutations(dObj, a[gene_samples,], b[gene_samples,], metric, cutoff, fun = perm_mode, blocs, num_perms)
				
				#Get p and q values of sampled genes
				curve_and_p <- get.filtered.genes(dObj, samples_perm, distance, gene_samples,metrics_and_cutoffs)		
				curve_genes <- curve_and_p[[1]]
				p_values_from_samples <- curve_and_p[[2]]
				q_values_from_samples<- curve_and_p[[3]]
			}
		}
	
		#Calculate p-values of good prospects 
		difconet.catf(dObj, "\nCalculating Permutations...\n")
		ptic = proc.time();
		permutations <-  get.permutations(dObj, a, b, metric, cutoff, fun = perm_mode, blocs, num_perms) 
		ptoc = proc.time();
		difconet.catf(dObj, "\nPermutations time:", (ptoc-ptic)["elapsed"], "sec\n");
		difconet.catf(dObj, "\nCalculating p-values...\n")
		prospects_p_values <- get.metric.p.values(dObj, distance, permutations, ifelse(speedup,curve_genes,rep(TRUE,nrow(distance))), metrics_and_cutoffs)
	
		#Estimate uncalculated p-values 
		if(speedup)
		{
			difconet.catf(dObj, "\nEstimating uncalculated p-values...\n")
			all_p_values <- estimate.metric.p.values(distance, gene_samples, curve_genes, p_values_from_samples, q_values_from_samples, prospects_p_values)
			est <- all_p_values[[1]] 
			p_vals <-  all_p_values[[2]]
			q_vals <- all_p_values[[3]]
			
			distance_p_values <- prepare.data(est, distance,p_vals,q_vals)			
		}
		else
		{
			distance_p_values <- cbind(rep(0, nrow(prospects_p_values)), prospects_p_values)
		}
		
		colnames(distance_p_values) <- c("Estimated", colnames(prospects_p_values))
		
		
		#Differential Expression p-values and q-values
		de.p <- apply(cbind(a,b), 1, function(x) { p<-NA;try(p<-wilcox.test(x[1:a_cols], x[(a_cols+1):both_cols])$p.value);p })
		de.q <- p.adjust(de.p, "fdr")
		info <- cbind(distance_p_values, de.p, de.q)
		colnames(info) <- c(colnames(distance_p_values),"expr.p", "expr.q")
		
		dObj$combstats[[index]] <- if (speedup) info else info[,-1]
		dObj$combdens[[index]] <- list(
			observed = apply(distance, 2, density, na.rm=TRUE),
			permutations = if (is.list(permutations) && length(permutations) > 0) lapply(permutations, function(x) density(x, na.rm=TRUE)) else list(density(permutations, na.rm=TRUE))
			)

		if (save_perm) {
			dObj$permutations[[index]] <- permutations
		}
		
	}
	dObj$global_p <- get.global.p.value(dObj)
	dObj$global_q <- p.adjust(dObj$global_p, method="fdr")
	
	return (dObj)
}

split.datset <- function(classes, all_data, labels_vector, comparisons)
{
	stages <- list()
	labels <- c()
	for (i in 1:length(classes)) {
		stages[[i]] <- data.matrix(all_data[,labels_vector  == classes[i]])
		labels <- c(labels, rep(classes[i],sum(labels_vector  == classes[i])))
		#difconet.catf(dObj, paste("Stage ", i ," :" , dim(stages[[i]])[2] ,sep = ""), "\n")
	}
	names(stages) <- classes

	combinations <- NA
	#Define the different stage pairwise experiments according to the @comparisons param
	if (is.character(comparisons) && length(comparisons) == 1 && comparisons %in% c("seq","all")) {
		if(comparisons == "seq"){
			#combinations <<- data.frame(lapply(strsplit(paste(labels_vector[1:(length(classes)-1)],"_",labels_vector[(2:length(classes))], sep = ""), "_"), as.numeric))
			combinations <- t(data.frame(classes[1:(length(classes)-1)], classes[(2:length(classes))], stringsAsFactors=FALSE))
		}else if(comparisons == "all"){
			#combinations <<- combn(1:length(classes), 2)
			combinations <- apply(combn(1:length(classes), 2), 2, function(x) classes[x])
		}
	}else{
		combinations <- data.frame(comparisons, stringsAsFactors=FALSE)
	}
	
	return (list(stages,combinations,labels))
}

get.global.p.value <- function(dObj) #data_final
{
	cs <- lapply(dObj$combstats, function(x) x[ ,grepl( "M[0-9\\.]+.p$" , colnames(x)), drop = FALSE])
	qs <- do.call(cbind, cs)
	if (ncol(qs) == 1) { global.p <- qs }
	else { global.p <- apply(qs, 1, function(x) pchisq( -2*sum(log(x), na.rm =  TRUE), df = 2*(length(x)-sum(is.na(x))), lower.tail=FALSE)) }

	return (global.p)
}

get.correlations.by.blocks <- function(dObj, a, b, metric, cutoff) {
 
	n <- nrow(a)
	total.genes <- ncol(a)
	distance <- list()
	
	#This is relevant to metric 5 only
	brks <- seq(-1, 1, length=201)
	brks[c(1,length(brks))] <- c(-1.1,1.1)
		
	if(1 %in% metric){
	
		difconet.catf(dObj, "M1: Sum(A>th)-Sum(B > th)...")
		tic = proc.time();
		a1 <- abs(a)
		b1 <- abs(b)
		
		for(ct in cutoff)
		{
			difconet.catf(dObj, "th=", ct, ",")
			#b.cnt <- (b1 > ct)
			#a.cnt <- (a1 > ct)
			#a.sum <- apply(a.cnt, 1, sum, na.rm=TRUE)
			#b.sum <- apply(b.cnt, 1, sum, na.rm=TRUE)
			a.sum <- rowSums(a1 > ct, na.rm=TRUE)
			b.sum <- rowSums(b1 > ct, na.rm=TRUE)
			distance <- c(distance,(data.matrix(abs(b.sum -a.sum))))
		}
		rm(a1)
		rm(b1)
		gc()
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");
	}
	if(2 %in% metric){
		difconet.catf(dObj, "M2: Kolmogorov-Smirnov D (A->B)...")
		tic = proc.time();
		a_cols <- ncol(a)
		b_cols <- ncol(b)
		both_cols <- a_cols + b_cols
		
		#distance <- c(distance, data.matrix(t(apply(cbind(a,b), 1, function(x) { unlist(ks.test(x[1:a_cols], x[(a_cols+1):both_cols])[c("statistic","p")])}))))
		#distance <- c(distance, data.matrix(t(sapply(1:nrow(a), function(x) { unlist(ks.test(a[x,], b[x,])[c("statistic","p")])}))))
		distance <- c(distance, data.matrix(sapply(1:nrow(a), function(x) { ks.test(a[x,], b[x,])$statistic })))
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");

	}
	if(3 %in% metric){
		difconet.catf(dObj, "M3: Sum(|A-B| > th)")
		tic = proc.time();
		ab <- abs(a - b)
		for(ct in cutoff)
		{
			difconet.catf(dObj, "th=", ct, ",")
			#distance <- c(distance,data.matrix(apply(ab, 1, function(x) sum(x > ct, na.rm=TRUE))))
			distance <- c(distance,data.matrix(rowSums(ab > ct, na.rm=TRUE)))
		}
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");
		rm(ab)
		gc()
	}
	if(4 %in% metric){
		difconet.catf(dObj, "M4: Euclidean Distance...")
		tic = proc.time();
        sq.diff <- (b-a)^2
		#distance <- c(distance,data.matrix(apply(sq.diff, 1, function(x) sqrt(sum(x, na.rm=TRUE))/total.genes)))
		distance <- c(distance,data.matrix(sqrt(rowSums(sq.diff, na.rm=TRUE))/total.genes))
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");
		rm(sq.diff)
		gc()
	}
	if(5 %in% metric){
		difconet.catf(dObj, "M5: Kullback Lieber...")
		tic = proc.time();
		M5 <- matrix(0, nrow=n, ncol=1)
		
		for (i in 1:n) {
			if (i %% 1000 == 0) difconet.catf(dObj, n-i, " left\n") else if (i %% 100 == 0) difconet.catf(dObj, ".")
		
			a.vector <- a[i,]
			b.vector <- b[i,]
			
			a.vector <- cut(a.vector, breaks=brks, labels=FALSE)
			b.vector <- cut(b.vector, breaks=brks, labels=FALSE)
			
			a.vector <- a.vector/sum(a.vector, na.rm=TRUE)
			b.vector <- b.vector/sum(b.vector, na.rm=TRUE)
			
			ok <- a.vector * b.vector 
			#We use 1/sum as an approximation to 0 
			ab <- ifelse(ok, log(a.vector/b.vector), log((1/sum(a.vector, na.rm=TRUE))/b.vector))
			ba <- ifelse(ok, log(b.vector/a.vector), log((1/sum(b.vector, na.rm=TRUE))/a.vector))
			M5[i] <- sum(a.vector * ab, na.rm=TRUE) + sum(b.vector * ba, na.rm=TRUE)
		}
		distance <- c(distance, M5)
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");
	}
	if(6 %in% metric)
	{
		difconet.catf(dObj, "M6: sqrt(.5|sign(b)*b^2-sign(a)*a^2|))^2.5...")
		tic = proc.time();
		bt <- 2.5
		a <- a*a*sign(a)
		b <- b*b*sign(b)
		dif <- sqrt(0.5*abs(b-a))^bt
		#distance <- c(distance, data.matrix(apply(dif, 1, sum)))
		distance <- c(distance, data.matrix(rowSums(dif, na.rm=TRUE)))
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");
		rm(dif)
	}
	if (8 %in% metric) {
		difconet.catf(dObj, "M8: (user func)...")
		tic = proc.time();
		distance <- c(distance, (data.matrix(dObj$metricfunc(dObj, a, b))))
		toc = proc.time();
		difconet.catf(dObj, (toc-tic)["elapsed"], "sec\n");
	}
	rm(a)
	rm(b)
	gc()
	distance <- do.call(cbind, distance)
	return (distance)	
}

difconet.get.cor.distance <- function(dObj, a, b, metric, cutoff, blocs)
{
	n <- nrow(a)
	metric_columns <- get.metrics.length(dObj, metric, cutoff)
	
	distance <- matrix(0, nrow=n, ncol=metric_columns)
   
	starts <- 1
	ends <-   min(blocs,n) #ifelse(blocs>n, n , blocs)
	check <- ceiling(n/blocs)
	#tic = proc.time();
	
	#a <- t(apply(a, 1, rank))
	#a = a - rowMeans(a)
	#a = a / sqrt(rowSums(a^2))
	
	#b <- t(apply(b, 1, rank))
	#b = b - rowMeans(b)
	#b = b / sqrt(rowSums(b^2))
	ta <- t(a)
	tb <- t(b)

	#toc = proc.time();
	#difconet.catf(dObj, "Preparation ", (toc-tic)["elapsed"], ", ");
   
	for(i in 1:check){    
		difconet.catf(dObj, paste( "Block " , i, ", ",starts,":",ends,"; "))
		tic = proc.time();
		#a.sub <- tcrossprod(a[starts:ends,], a)
		#b.sub <- tcrossprod(b[starts:ends,], b)
		a.sub <- dObj$corfunc(ta[,starts:ends], ta)
		b.sub <- dObj$corfunc(tb[,starts:ends], tb)
		toc = proc.time();
		difconet.catf(dObj, "Corr", (toc-tic)["elapsed"], "sec\n");
		distance[starts:ends,] <- get.correlations.by.blocks(dObj, a.sub, b.sub, metric, cutoff)
		starts <- starts+blocs 
		ends <- ifelse((ends+blocs)>n, n , ends+blocs ) 
	}
	colnames(distance) <- paste("M", get.metric.cutoffs(dObj, metric, cutoff), sep = "")
	return(distance)
 }

get.metrics.length <- function(dObj, metric, cutoff)
{
	cutoff.count <- (1 %in% metric) + (3 %in% metric)
	len <- length(metric) + (cutoff.count)*length(cutoff) - cutoff.count
}

get.metric.cutoffs <- function(dObj, metric, cutoff)
{
	if((1 %in% metric) || (3 %in% metric))
	{
		pos <- c(which(metric == 1), which(metric == 3))
		vals <- list()
		i <- 1
		for(m in pos)
		{
			vals[[i]] <-  (metric[m] + cutoff)
			i <- i+1
		}
		m.names <- insert.at(as.character(metric), pos[1], vals) 
		pos <- c(which(m.names == 1), which(m.names == 3))
		m.names <- m.names[-pos]	
	}
	else
		m.names <- metric
	return(m.names)
}

get.permutations <- function(dObj, a, b, metric, cutoff, fun = "columns", blocs, num_perm)
{
	a_cols <- ncol(a)
	b_cols <- ncol(b)
	both_cols <- a_cols + b_cols
	metrics_and_cutoffs <- get.metric.cutoffs(dObj, metric, cutoff)
	
	permutations <-  list()
	nrows <- nrow(a) #ifelse(length(selection)==0, nrow(a), length(selection))
	for(ms in 1:length(metrics_and_cutoffs))
	{
		permutations[[ms]] <- matrix(0, nrow=nrows, ncol=num_perm)
	}
	names(permutations) <- paste0("M",metrics_and_cutoffs)

	for(i in 1:num_perm){
		difconet.catf(dObj, "--Iteration ", i, "...")
		tic = proc.time();
		if (fun == "columns") {
			#permutaciones de etiquetas
			blob <- cbind(a, b)
			x <- sample(1:both_cols, both_cols, replace=F)
			blob <- blob[,x]
			a_new <- blob[,1:a_cols]
			b_new <- blob[,(a_cols+1):both_cols]
		} else if (fun == "rows") {
			#Permutaciones de etiquetas por gen
			blob <- cbind(a, b)
			xt <- t(apply(blob, 1, sample))
			a_new <- xt[,1:a_cols]
			b_new <- xt[,(a_cols+1):both_cols]
		} else if (fun == "rows.class") {
			#Permutations de etiquetas por gen en cada clase
			a_new <- t(apply(a, 1, sample))
			b_new <- t(apply(b, 1, sample))
		} else if (fun == "all") {
			blob <- cbind(a, b)
			xt <- matrix(sample(blob), ncol=ncol(blob))
			a_new <- xt[,1:a_cols]
			b_new <- xt[,(a_cols+1):both_cols]
		}

		new.distance <- difconet.get.cor.distance(dObj, a_new, b_new, metric, cutoff, blocs)
		for(ms in 1:length(metrics_and_cutoffs))
		{
			permutations[[ms]][,i] <- new.distance[,ms]
		}
		rm(a_new)
		rm(b_new)
		rm(blob)
		gc()
		toc = proc.time();
		difconet.catf(dObj, "-- ", i, (toc-tic)["elapsed"], "sec\n");
	}	

	return (permutations)
}

get.representative.genes <- function(dObj, dist_vector, quant)
{
	genes_by_quantile <- list()
	sample_by_quantile <- ceiling(length(dist_vector)*0.01)
	for(i in 1:10)
	{
		lower_quantile <- quant[i] 
		upper_quantile <- quant[i+1]
		quantile_filter <- dist_vector > lower_quantile & dist_vector < upper_quantile
		genes_by_quantile[[i]] <- sample(which(quantile_filter), sample_by_quantile)	
	}
	genes_by_quantile <- do.call(c, genes_by_quantile)
	return(genes_by_quantile)
}

get.filtered.genes <- function(dObj, samples_perm, distance, gene_samples,metric)
{
	p_values_from_samples <- list()
	q_values_from_samples<- list()
	curve_genes <- list()
	for(m in 1:length(samples_perm))
	{
		calculated <- get.distance.p.values(dObj, distance[gene_samples,m], samples_perm[[m]])
		p_values_from_samples[[m]] <- calculated[,2]
		q_values_from_samples[[m]] <- calculated[,3]
		interest_point <- min((distance[gene_samples,m])[p_values_from_samples[[m]] <= (1/dObj$num_perms)])
		difconet.catf(dObj, paste("Interest Point @ M", metric[m], " : ", interest_point, "\n", sep = ""))
		curve_genes[[m]] <- distance[,m] >= interest_point
		difconet.catf(dObj, paste("Number of genes @ M", metric[m], " : ", sum(curve_genes[[m]]), "\n", sep = ""))
	}
	
	curve_genes <- do.call(cbind, curve_genes)
	curve_genes <- apply(curve_genes, 1, function(x) any(x))
	p_values_from_samples <- do.call(cbind, p_values_from_samples)
	q_values_from_samples<- do.call(cbind, q_values_from_samples)
	return(list(curve_genes,p_values_from_samples, q_values_from_samples))
}

get.metric.p.values <- function(dObj, distance, permutations, curve_genes, metric)
{
	distance_p_values <- list()

	for(m in 1:ncol(distance))
	{
		difconet.catf(dObj, "p-values metric: " , metric[m], "\n")
		distance_p_values[[m]] <- get.distance.p.values(dObj, distance[curve_genes,m],permutations[[m]])
		colnames(distance_p_values[[m]]) <- paste("M", metric[m],c(".dist",".p", ".q"), sep ="")
	}
	distance_p_values <- do.call(cbind, distance_p_values)
	return(distance_p_values)
}


estimate.metric.p.values <-function(distance, gene_samples, curve_genes, p_values_from_samples, q_values_from_samples, distance_p_values)
{
	available_genes <- unique(c(gene_samples, which(curve_genes)))
	estimated_p_values  <- list()	
	estimated_q_values  <- list()
	for(m in 1:ncol(p_values_from_samples)){
		lm_distance <- distance[available_genes, m]
		
		#Bind both sample and goos prospect genes' p-values and their correction
		p_vals <- rep(0, dim(distance)[1])
		p_vals[gene_samples] <- p_values_from_samples[,m]
		p_vals[which(curve_genes)] <- distance_p_values[,(3*m)-1]
		lm_p_vals <- p_vals[available_genes]
		
		
		q_vals <- rep(0, dim(distance)[1])
		q_vals[gene_samples] <- q_values_from_samples[,m]
		q_vals[which(curve_genes)] <- distance_p_values[,(3*m)]
		lm_q_vals <- q_vals[available_genes]
		
		exp_model <-lm(log(lm_p_vals)  ~ lm_distance)
		remaining <- distance[-available_genes, m]
		p_vals[-available_genes] <- exp(predict(exp_model,list(lm_distance=remaining)))
		p_vals[p_vals > 1.00] <- 1.0
		estimated_p_values[[m]] <- p_vals
		
		exp_model <-lm(log(lm_q_vals)  ~ lm_distance)
		q_vals[-available_genes] <- exp(predict(exp_model,list(lm_distance=remaining)))
		q_vals[q_vals > 1.00] <- 1.0
		estimated_q_values[[m]] <- q_vals
	}
	estimated_p_values <- do.call(cbind, estimated_p_values)
	estimated_q_values <- do.call(cbind, estimated_q_values)
	estimated <- rep(FALSE, dim(distance)[1])
	estimated[-available_genes] <- TRUE
	
	return(list(estimated, estimated_p_values, estimated_q_values))
}

#Bind different data into a single dataframe
prepare.data <- function(est, distance,p_vals,q_vals)
{
	all_info <-  data.frame(est)
	for(m in 1:ncol(distance))
	{
		all_info <- cbind(all_info,distance[,m], p_vals[,m], q_vals[,m])
	}
	return(all_info)
}

difconet.plot.gene.correlations <- function(dObj, gene, stages=1:length(dObj$stages.data), type=c("density","scatter")[1], main=rownames(dObj$stages.data[[1]])[gene], legends=TRUE, plot=TRUE, ...)
{
	corSt <- NULL
	for (i in stages) {
		a <- dObj$stages.data[[i]]
		xcor <- dObj$corfunc(t(a), t(a[gene, , drop = FALSE])) #, method = "spearman"
		corSt <- cbind(corSt, xcor)
	}
	colnames(corSt) <- names(dObj$stages.data)[stages]
	type <- match.arg(type, c("density","scatter"))
	if (type == "density") {
		if (plot) {
			plot(density(corSt, na.rm=TRUE), type="n", xlab="", ylab="Density", main=main, ...)
			for (i in 1:length(stages)) {
				lines(density(corSt[,i], na.rm=TRUE), col=i)
			}
			if (legends) legend("topright", dObj$labels[stages], col=stages, lty=1, bg=0, box.col=0)
		}
	} else {
		#require(grDevices)
		rbPal <-   colorRampPalette(c("#DDDDDD", "#00008F", "#0000DD", "#007FDD", "#00DDDD",
				"#7FDD7F", "#DDDD00", "#DD7F00", "#DD0000", "#DD30DD"))
		rbPal <- rbPal(1000)
		rbPal2 <- paste0(rbPal, rep(c(1,3,5,7,9,"B","C","D","E","F"), each=100), "0")
		panelito <- function(a, b, ...) {
			xdif <- (a-b)^2
			points(a, b, 
				xlim = c(-1,1), ylim = c(-1,1),
				col=rbPal[cut(xdif,breaks=c(seq(0,1,length.out=1000),Inf),labels=FALSE)], 
				pch=20)
			abline(0,1, lty = 2)			
		}
		panelote <- function(a, b, ...) {
			xdif <- (a-b)^2
			points(a, b, 
				xlim = c(-1,1), ylim = c(-1,1),
				col=rbPal2[cut(xdif,breaks=c(seq(0,1,length.out=1000),Inf),labels=FALSE)], 
				pch=20)
			abline(0,1, lty = 2)			
		}
		if (plot) {
			if (length(stages) > 2) {
				pairs(corSt, gap=0, 
					xlim = c(-1,1), ylim = c(-1,1),
					upper.panel=panelote,
					lower.panel=panelito,
					pch=20, main=main)
				#xlab = dObj$labels[stages[1]], ylab = dObj$labels[stages[2]], main=main)
			} else {
				xdif <- (corSt[,1]-corSt[,2])^2
				plot(corSt[,1], corSt[,2], 
					xlim = c(-1,1), ylim = c(-1,1),
					col=rbPal[cut(xdif,breaks=c(seq(0,1,length.out=1000),Inf),labels=FALSE)], 
					pch=20, xlab = dObj$labels[stages[1]], ylab = dObj$labels[stages[2]], main=main)
				abline(0,1, lty = 2)
			}
		}
	}
	invisible(corSt)
}

get.distance.p.values <- function(dObj, distance, permutations)
{
	p_values <- matrix(0, nrow=length(distance), ncol=2)
	
	if (! dObj$use_all_perm) {
		P <- ncol(permutations)
		for (i in 1: length(distance)) {
			NP <- sum(permutations[i,] > distance[i])
			p_values[i,1] = (NP + 1)/P
		}
	} else {
		P <- length(permutations)
		for (i in 1: length(distance)){
			NP <- sum(permutations > distance[i])
			p_values[i,1] = (NP + 1)/P
		}
	}
	
	p_values[,2] = p.adjust(p_values[,1], "fdr")
	
	distance_p_values  <- cbind(distance, p_values)
	
	return (distance_p_values)
}

difconet.plot.histograms.heatmap2 <- function(dObj, genes=1:10, stages=1:length(dObj$stages.data), qprobs=c(0,.50,.975,.995), ...)
{
	corSt <- list()
	csc <- c()
	for (i in stages) {
		a <- dObj$stages.data[[i]]
		xcor <- t(dObj$corfunc(t(a), t(a[genes, , drop = FALSE])))
		corSt[[i]] <- t(get.histogram(dObj, xcor))
		csc <- c(csc, rep(i, ncol(corSt[[i]])))
	}
	names(corSt) <- names(dObj$stages.data)[stages]
	corH <- do.call(cbind, corSt)
	q <- quantile(corH, probs = qprobs)
	breakz = c(seq(q[1],q[2],length.out=12),
		seq(q[2],q[3],length.out=24)[-1],
		seq(q[3],q[4],length.out=12)[-1])
	Colors = c(	'#A6A6A6','#FFFFFF','#CCCCE8','#9999D1','#6666B9','#3333A2',
				'#00008B','#00146F','#002853','#003C38','#00501C','#006400',
				'#006400','#187300','#308100','#489000','#619F00','#79AE00',
				'#91BC00','#AACB00','#C2DA00','#DAE900','#F2F700','#FAF300',
				'#EFDB00','#E4C300','#D8AA00','#CD9200','#C27A00','#B86200',
				'#AC4900','#A13100','#961900','#8B0000','#8B0000','#A20000',
				'#B90000','#D10000','#E80000','#FF0000','#EC0630','#D90D60',
				'#C61390','#B31AC0','#A020F0')
	#require(gplots)
	heatmap.2(corH, Colv=FALSE, breaks=breakz, col=Colors, colsep=0:(ncol(corH)/50)*50, margins=c(5,10), 
		sepcolor="black", trace="none", ColSideColors=as.character(csc), dendrogram="row", ...)
	title("Correlation Density")
	cat("Color Order:",paste(dObj$labels[stages], collapse="; "),"\n")
}

get.histogram <- function(dObj, dataset, brks =seq(-1,1,by=0.01)) #a-b
{
	histograms <- apply(dataset, 1, function(x) hist(x, breaks=brks, plot=FALSE)$counts)
	#histograms <- matrix(0, nrow=nrow(dataset), ncol=200)
	#for(i in 1:nrow(dataset)){
	#	histograms[i,] <- hist(dataset[i,],  breaks= brks, plot = FALSE)$counts
	#}
	return  (histograms)
}

insert.at <- function(a, pos, ...){
    result <- vector("list",2*length(pos)+1)
    result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
    result[c(FALSE,TRUE)] <- list(...)
    return(unlist(result))
}


