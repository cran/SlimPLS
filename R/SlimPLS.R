library(e1071)
library(randomForest)
library(bioDist)
library(zipfR)

setClass("expMat", representation(data = "matrix", 
								  features = "vector",
				  				  samples = "vector",
								  classes = "vector" ) )
				  
setClass("featureSet", representation(features = "vector", t_mat_learn = "matrix", w_mat = "matrix", 
									  p_mat = "matrix", q_mat = "matrix", is_pls = "logical" ) )

setOldClass(c("svm","randomForest"))

setClass("svmOrRF")

setIs("svm", "svmOrRF")
setIs("randomForest", "svmOrRF")

setClass("classificationModel", representation( model_type="character", model = "svmOrRF", knn_k = "numeric",
												features = "featureSet", pls_use_components = "logical") ) 

		
##
# Read an expression matrix from a file into a expMat obeject.
# Feature labels are expected to be in first column
# Sample labels are expected to be in first row
# Optionally, the matrix can also hold classes for the different samples. If exist, they
# are expected to be in the second row
# File should be tab delimited.
##
readExpMat <- function(file_name="", with_classes = TRUE)
{
	mat <- as.matrix(read.table(file_name,sep="\t")) 

	end_row <- nrow(mat)
	end_col <- ncol(mat)
	
	if (with_classes == TRUE ) {
		start_row <- 3
		classes_vector <- as.vector(mat[2,2:end_col])
	}
	else
	{
		start_row <- 2
		classes_vector <- c(0)
	}
		
	
	sub_mat <- mat[start_row:end_row,2:end_col]
	sub_mat <- matrix(as.numeric(sub_mat),nrow(sub_mat),ncol(sub_mat))
	sub_mat <- t(sub_mat)
	
	samples_vector <- as.vector(mat[1,2:end_col])
	features_vector <- as.vector(mat[start_row:end_row,1])
	
	rownames( sub_mat ) <- samples_vector
	colnames( sub_mat ) <- features_vector
			
	exp_mat<-new("expMat", data=sub_mat, features = features_vector, samples = samples_vector, classes = classes_vector )
	
	return(exp_mat)	
}

##
# Selects the best features using the given select_method from the given expression matrix exp_mat
# The expression matrix is expected to be ordered with all the samples belonging to the first class
# grouped at the beginning, and all the samples belonging to the second class grouped together
# following them. The number of samples from the two classes are given in num_class_a and
# num_class_b
# select_method can get the following values:
# CORREL - Pearson correlation score
# TTEST - t test score with p-value threshold,
# FC - fold change score
# GOLUB - Golub criterion score
# MI - mutual information score
# The function selects the top scoring features according to the given score, the maximal number
# of features to be selected is given in num_features. 
##
selectFeatures <- function(select_method, exp_mat, num_class_a = 0 , num_class_b = 0, 
		class_a="", class_b="", num_features )
{
	if( num_class_a > 0 )
	{
		train_mat = exp_mat@data[1:(num_class_a+num_class_b),]
	}
	else
	{
		class_a_indices = which(exp_mat@classes==class_a)
		class_b_indices = which(exp_mat@classes==class_b)
		all_classes = c(class_a_indices, class_b_indices)
		num_class_a = length(class_a_indices)
		num_class_b = length(class_b_indices)
		train_mat = exp_mat@data[all_classes,]
	}
	
	filter_type <- PLS_dispatch_filter_type( select_method )
	if (select_method == 'TTEST')
	{
		dec = FALSE
	}
	else
	{
		dec = TRUE
	}
	index_vec <- get_top_features_using_filter(train_mat, num_class_a, num_class_b,
											   num_features, filter_type, decreasing = dec)
	
	features_vector <- (exp_mat@features)[ index_vec ]
	
	feature_set <- new("featureSet", features = features_vector, is_pls = FALSE )
	return(feature_set)
}

##
# Selects the best features using SlimPLS algorithm from the given expression matrix exp_mat
# The expression matrix is expected to be ordered with all the samples belonging to the first class
# grouped at the beginning, and all the samples belonging to the second class grouped together
# following them. The number of samples from the two classes are given in num_class_a and
# num_class_b
# The function selects up to num_features features, grouped together is components of the given
# component_size. If p_value_threshold is not 0, it selects only features with a p-value lower than
# the given threshold.
# The features can later be used as components (linear combination of features) or simply as a list 
# of features.
##
selectFeaturesSlimPLS <- function(exp_mat, num_class_a=0, num_class_b=0, 
		class_a="", class_b="", num_features = 50, component_size = 25, p_value_threshold = 0)
{
	
	if( num_class_a > 0 )
	{
		train_mat = exp_mat@data[1:(num_class_a+num_class_b),]
	}
	else
	{
		class_a_indices = which(exp_mat@classes==class_a)
		class_b_indices = which(exp_mat@classes==class_b)
		all_classes = c(class_a_indices, class_b_indices)
		num_class_a = length(class_a_indices)
		num_class_b = length(class_b_indices)
		train_mat = exp_mat@data[all_classes,]
	}
	y_res <- c(rep(1,num_class_a),rep(2,num_class_b))
	
	if( p_value_threshold == 0 )
	{
		num_comp <- PLS_find_num_components(train_mat, y_res, num_features, component_size)
	}
	else
	{
		num_comp <- PLS_find_num_components_using_p_val(train_mat, y_res, num_features, p_value_threshold)
	}
	index_vec <- PLS_find_last_index_vector(num_comp)

	t_mat_learn_cached <- PLS_core_read_file("PLS_t_mat_learn.txt")
	w_mat_cached <- PLS_core_read_file("PLS_w_mat.txt")
	q_mat_cached <- PLS_core_read_file("PLS_q_mat.txt")
	p_mat_cached <- PLS_core_read_file("PLS_p_mat.txt") 
	
	features_vector <- (exp_mat@features)[ c(index_vec==1) ]
	w_mat_cached <- w_mat_cached[c(index_vec==1),]
	
	rownames( w_mat_cached ) <- features_vector
	names( features_vector ) <- features_vector
	feature_set <- new("featureSet", features = features_vector, t_mat_learn = t_mat_learn_cached, w_mat = w_mat_cached,
									 q_mat = q_mat_cached, p_mat = p_mat_cached, is_pls = TRUE)
	
	return(feature_set)
	
}

trainClassifier <- function(learn_method, exp_mat, num_class_a=0, num_class_b=0,
		class_a="", class_b="", feature_set, use_components_as_features = FALSE )
{
	
	if( use_components_as_features == TRUE)
	{
		PLS_method = PLS_get_t_factors
	}
	else
	{
		PLS_method = PLS_get_features_by_w
	}
	
	if( num_class_a > 0 )
	{
		test_mat = exp_mat@data[1:(num_class_a+num_class_b),]
	}
	else
	{
		class_a_indices = which(exp_mat@classes==class_a)
		class_b_indices = which(exp_mat@classes==class_b)
		all_classes = c(class_a_indices, class_b_indices)
		num_class_a = length(class_a_indices)
		num_class_b = length(class_b_indices)
		test_mat = exp_mat@data[all_classes,]
	}

	y_res <- c(rep(1,num_class_a),rep(2,num_class_b))
	
	# if this is a SlimPLS selected features
	if ( feature_set@is_pls )
	{
		num_features <- length(exp_mat@features)
		num_components <- ncol(feature_set@w_mat)
		w_mat <- matrix( rep(0, num_components*num_features), ncol = num_components, nrow = num_features ) 
		for(i in 1:num_features)
		{
			
			if( !is.na(feature_set@features[exp_mat@features[i]]))
			{
				w_mat[i,] <- feature_set@w_mat[ exp_mat@features[i],]
				
			}
		}
		PLS_core_write_file(w_mat, file = "PLS_w_mat.txt")
		x_mat <- as.matrix(test_mat, nrow= num_class_a+num_class_b)
		new_learn_x_mat <- PLS_method(x_mat, learn=TRUE, num_c=num_components)
	}
	else
	{
		new_learn_x_mat = as.matrix(test_mat, nrow= num_class_a+num_class_b)
		# todo: project matrix on filtered feature set
	}
	
	
	
	if(learn_method == "SVM_LINEAR")
	{
		c_par_linear <- find_c_regularization(new_learn_x_mat, factor(y_res), kernel="linear", type="C-classification")
		res_learn <- svm(new_learn_x_mat, factor(y_res), kernel="linear", type="C-classification", cost=c_par_linear)
		cls_model <- new("classificationModel", model_type = learn_method, model = res_learn, features = feature_set,
				pls_use_components = use_components_as_features )
	}
	if(learn_method == "SVM_RADIAL")
	{
		c_par_radial <- find_c_regularization(new_learn_x_mat, factor(y_res), kernel="radial", type="C-classification")
		res_learn <- svm(new_learn_x_mat, factor(y_res), kernel="radial", type="C-classification", cost=c_par_radial)
		cls_model <- new("classificationModel", model_type = learn_method, model = res_learn, features = feature_set,
				pls_use_components = use_components_as_features)
	}
	if(learn_method == "RANDOM_FOREST")
	{
		res_learn <- randomForest(new_learn_x_mat, factor(y_res), ntree=1500)
		cls_model <- new("classificationModel", model_type = learn_method, model = res_learn, features = feature_set,
				pls_use_components = use_components_as_features)
	}
	
	return (cls_model)
	
}

getClassification <- function(exp_mat, model)
{
	if( model@pls_use_components == TRUE)
	{
		PLS_method = PLS_get_t_factors
	}
	else
	{
		PLS_method = PLS_get_features_by_w
	}
	learn_method <- model@model_type
	feature_set <- model@features
	
	x_mat <- exp_mat@data
	
	
	# if this is a SlimPLS selected features
	if ( feature_set@is_pls )
	{
		num_features <- length(exp_mat@features)
		num_components <- ncol(feature_set@w_mat)
		w_mat <- matrix( rep(0, num_components*num_features), ncol = num_components, nrow = num_features ) 
		for(i in 1:num_features)
		{
			
			if( !is.na(feature_set@features[exp_mat@features[i]]))
			{
				w_mat[i,] <- feature_set@w_mat[ exp_mat@features[i],]
				
			}
		}
		PLS_core_write_file(w_mat, file = "PLS_w_mat.txt")
		
	}
	
	n_samples <- nrow(exp_mat@data)
#	print(n_samples)
	res_vec <- rep(0,n_samples)
		
	if(learn_method == "SVM_LINEAR" || learn_method == "SVM_RADIAL" || learn_method == "RANDOM_FOREST")
	{
		
		for(i in 1: n_samples)
		{
			if ( feature_set@is_pls )
			{
				
				PLS_core_predict(x_mat[i,], num_components, wf = TRUE)
				
				new_predict_x_mat <- PLS_method(x_mat[i,], learn=FALSE, num_c=num_components)
			}
			else
			{
				new_predict_x_mat <- x_mat[i,]
			}
			res_pred <- predict(model@model, t(as.matrix(new_predict_x_mat)))
			res_vec[i] <- res_pred
		}
	}
	names(res_vec)=exp_mat@samples
	return(res_vec)
}

##
# Scoring functions for filter based selection method. Each function scores a feature
# according to its value over the different samples, and the labels given to each of the
# samples.
# function prototype is function(mat, j, size_w1, size_w2, ...)
# mat is the expression matrix. j is the number of feature to test. size_w1 and size_w2
# are the number of samples belonging to class 1 and 2.
##

mutual_information_score <- function(mat, j, size_w1, size_w2, ...)
{
	c_vector <- c( rep(1,size_w1), rep(2,size_w2) )
	f_vector <- mat[,j]
	ref_mat <- rbind(c_vector,f_vector)
	res <- as.double(mutualInfo(ref_mat))
	return(res)
}

golub_criterion_score <- function(mat, j, size_w1, size_w2, ...)
{
	vector <- mat[,j]
	g_w1 <- vector[1:size_w1]
	g_w2 <- vector[(size_w1+1):(size_w1+size_w2)]
	tmp_numerator <- mean(g_w1) - mean(g_w2)
	tmp_denominator <- sd(g_w1) + sd(g_w2)
	res <- abs(tmp_numerator/tmp_denominator)
	return(res)
}

t_test_p_val_score <- function(mat, j, size_w1, size_w2, get_pval = TRUE, ...)
{
	vector <- mat[,j]
	if (get_pval == TRUE) {
		res <- (t.test(vector[1:size_w1],vector[(size_w1+1):(size_w1+size_w2)]))$p.value	
	}
	else {
		res <- (t.test(vector[1:size_w1],vector[(size_w1+1):(size_w1+size_w2)]))$statistic
	}
	return(res)
}

PLS_calc_p_val_from_cor <- function(r, n_samples)
{
	df <- n_samples - 2
	tstat <- (r * sqrt(df))/(sqrt(1-(r^2)))
	#print(tstat)
	return(PLS_calc_p_val_from_ttest(tstat,df))	
}

PLS_calc_p_val_from_ttest <- function(tstat, df)
{
	res <- Rbeta(df/(df+(tstat^2)),0.5*df,0.5,lower=TRUE,log=FALSE)
	#or
	#x <- (tstat + sqrt(tstat^2+df))/(2*sqrt(tstat^2 + df))
	#res <- Rbeta(x, df/2, df/2, lower=TRUE, log=FALSE)
	#res <- 2*(1-res)
	return(res)
}

max_expression_level <- function(mat, j, size_w1, size_w2, ...)
{
	vector <- mat[,j]
	res <- max(vector)
	return(res)
}



fold_change_score <- function(mat, j, size_w1, size_w2, norm=TRUE, ...)
{
	vector <- mat[,j]
	res <- mean(vector[1:size_w1]) / mean(vector[(size_w1+1):(size_w1+size_w2)]) 
	if (norm == TRUE) {
		if (res < 1) {
			res <- 1/res
		}
	}
	return(res)
}

correlation_score <- function(mat, j, size_w1, size_w2, use_abs=TRUE, ...)
{
	class_vector <- c( rep(1,size_w1), rep(2,size_w2) )
	cor_res <- cor(mat[,j], class_vector)
	if (use_abs == TRUE) {
		cor_res <- abs(cor_res)
	}
	return(cor_res)
}


##
# returns the appropriate scoring function according to the given string
# filter_str can get the following values:
# CORREL - Pearson correlation score
# TTEST - t test score with p-value threshold,
# FC - fold change score
# GOLUB - Golub criterion score
# MI - mutual information score
##
PLS_dispatch_filter_type <- function(filter_str)
{
	ret_type <- switch(filter_str,
			CORREL = correlation_score,
			TTEST = t_test_p_val_score,
			FC = fold_change_score,
			GOLUB = golub_criterion_score,
			MI = mutual_information_score)
	if (length(ret_type) == 0) {
		ret_type <- max_expression_level
	}
	return(ret_type)
}

##
# returns the top f_num features, using the given scoring function test_alg
##
get_top_features_using_filter <- function(org_mat, size_w1, size_w2, f_num, test_alg, decreasing=FALSE,...)
{
	score_vector <- get_score_vector(org_mat, test_alg, size_w1, size_w2, ...)
	index_list <- PLS_get_index_vector(score_vector, f_num, decreasing)
	return(index_list) 
}

##
# get vector of scores, one score for each feature, using the given scoring function
# test_alg
##
get_score_vector <- function(mat, test_alg, size_w1, size_w2, ...)
{
	score_vector <- c(0, ncol(mat))
	
	for (j in 1:ncol(mat)) {
		score_vector[j] <- test_alg(mat, j, size_w1, size_w2, ...)
	}
	return(score_vector)
}

##
#
# finding top k (top min or top max)
#
##
PLS_get_index_vector <- function(vec, num_f, decreasing=FALSE, ...)
{
	num_found = 0
	l_vec <- length(vec)
	ret_vec <- rep(FALSE, l_vec)
	
	if (num_f > l_vec) {
		return(!ret_vec)
	}
	
	if (decreasing == TRUE) {
		min_or_max <- max
		n_value <- -Inf
	}
	else {
		min_or_max <- min
		n_value <- Inf
	}
	
	while (num_found < num_f) {
		cur_val <- min_or_max(vec)
		cur_index <- (vec == cur_val)
		vec[cur_index] <- n_value
		ret_vec <- ret_vec | cur_index
		num_found <- length(ret_vec[ret_vec == TRUE]) 
	}
	return(ret_vec)
}

##
# PLS functions
##

PLS_find_num_components <- function(x_mat, y_res, max_num, k_one, multi=2)
{
	f_vector <- rep(FALSE, ncol(x_mat))
	total_num <- 0
	num_iter <- as.integer((max_num/k_one)*multi)
	features_selected <- PLS_core_learning(x_mat, y_res, num_iter, w_adjustment=PLS_adj_top_k_features, k_features=k_one)
	w_mat <- PLS_core_read_file("PLS_w_mat.txt")
	
	for (i in 1:(ncol(w_mat))) {
		tmp_vec <- w_mat[,i]
		tmp_vec[tmp_vec==0] <- FALSE
		tmp_vec[tmp_vec!=0] <- TRUE
		f_vector <- f_vector | tmp_vec
		if ((length(f_vector[f_vector==TRUE]) == max_num) || ((i == 1)&&(length(f_vector[f_vector==TRUE]) > max_num))){
#			print("Exact match")
			PLS_core_write_file(w_mat[,1:i], file = "PLS_w_mat.txt")
			return(i)
		}
		else if (length(f_vector[f_vector==TRUE]) > max_num) {
#			print("NOT Exact match")
			#tvec <- prev_f_vector & tmp_vec
			#tmp_vec[t_vec] <- FALSE
			#tmp_w <- (w_mat[,1])[t_vec]
			#needed_f <- max_num - length(prev_f_vector[prev_f_vector==TRUE])
			#tmp_w <- PLS_adj_top_k_features(tmp_w, needed_f)
			#w_mat[,i] <- tmp_w
			#PLS_core_write_file(w_mat[,1:i], file = "PLS_w_mat.txt")
			#return(i)
			PLS_core_write_file(w_mat[,1:(i-1)], file = "PLS_w_mat.txt")
			return(i-1)
		}
		#prev_f_vector <- f_vector
	} 
	return(PLS_find_num_components(x_mat,y_res,max_num,k_one,multi+1))
}

PLS_find_num_components_using_p_val <- function(x_mat, y_res, max_num, th)
{
	res <- PLS_core_k_learning(x_mat, y_res, max_num, th)
	PLS_core_write_file(res, "PLS_num_f.txt")
	w_mat <- PLS_core_read_file("PLS_w_mat.txt")
	return(ncol(w_mat))
}

PLS_core_learning <- function(x_matrix, y_response, num_iterations, w_adjustment = NULL, fs = FALSE, wf = TRUE, get_cor = FALSE, file_name="", ...)
{
	features_selected <- rep(0, ncol(x_matrix))
	q_mat <- rep(0, num_iterations)
	p_mat <- matrix(0, ncol(x_matrix), num_iterations)
	w_mat <- matrix(0, ncol(x_matrix), num_iterations)
	t_mat <- matrix(0, nrow(x_matrix), num_iterations)
	cor_vec <- rep(0, num_iterations)
	
	if (length(w_adjustment) == 0) {
		w_adjustment <- PLS_core_default_w_adjustment 
	}
	x_matrix <- PLS_core_scale_matrix(x_matrix)
	y_response <- PLS_core_scale_vector(y_response)
	y_bul <- rep(0,length(y_response)) #
	y_target <- y_response
	n_samples <- nrow(x_matrix)	
	
	for (i in 1:num_iterations) {
		
		w_loading <- PLS_core_get_w_loading(x_matrix, y_response)
		w_loading <- w_adjustment(w_loading, x_matrix=x_matrix, y_response=y_response, ...)
		w_mat[,i] <- w_loading
		
		if (fs == TRUE) {
			features_selected <- PLS_core_update_features_selected(w_loading, features_selected)
		}		
		
		t_factor <- x_matrix %*% w_loading
		t_mat[,i] <- t_factor
		p_loading <- PLS_core_get_loadings(x_matrix, t_factor)
		p_mat[,i] <- p_loading
		q_loading <- PLS_core_get_loadings(y_response, t_factor)
		q_mat[i] <- q_loading
		e_error <- x_matrix - t_factor %*% t(p_loading)
		y_est <- t_factor %*% t(q_loading)
		f_error <- y_response - y_est
		
		if (get_cor == TRUE) {
			cc <- as.numeric(cor(t_factor, y_target))
			cor_res <- PLS_calc_p_val_from_cor(cc,n_samples)
			cor_vec[i] <- cor_res
		}
		
		x_matrix <- e_error
		y_response <- f_error
	}
	
	if (wf == TRUE) {
		PLS_core_write_file(q_mat, file = "PLS_q_mat.txt")
		PLS_core_write_file(p_mat, file = "PLS_p_mat.txt")
		PLS_core_write_file(w_mat, file = "PLS_w_mat.txt")
		PLS_core_write_file(t_mat, file = paste("PLS_t_mat_learn",file_name,".txt",sep=""))		
	}
	
	if (get_cor == TRUE) {
		return(cor_vec)
	}
	return(features_selected)
}

PLS_core_k_dist_learning <- function(x_matrix, y_response, w_dist, fs=FALSE, wf=TRUE, file_name="")
{
	n_samples <- nrow(x_matrix)
	n_features <- ncol(x_matrix)
	num_iterations = length(w_dist)
	features_selected <- rep(0, n_features)
	q_mat <- rep(0, num_iterations)
	p_mat <- matrix(0, n_features, num_iterations)
	w_mat <- matrix(0, n_features, num_iterations)
	t_mat <- matrix(0, n_samples, num_iterations)
	
	x_matrix <- PLS_core_scale_matrix(x_matrix)
	y_response <- PLS_core_scale_vector(y_response)
	y_target <- y_response
	
	for (i in 1:num_iterations) {
		
		w_loading <- PLS_core_get_w_loading(x_matrix, y_response)
		w_loading <- PLS_adj_top_k_features(w_loading, k_features=w_dist[i],x_matrix=x_matrix, y_response=y_response)
		w_mat[,i] <- w_loading
		
		if (fs == TRUE) {
			features_selected <- PLS_core_update_features_selected(w_loading, features_selected)
		}		
		
		t_factor <- x_matrix %*% w_loading
		t_mat[,i] <- t_factor
		p_loading <- PLS_core_get_loadings(x_matrix, t_factor)
		p_mat[,i] <- p_loading
		q_loading <- PLS_core_get_loadings(y_response, t_factor)
		q_mat[i] <- q_loading
		e_error <- x_matrix - t_factor %*% t(p_loading)
		y_est <- t_factor %*% t(q_loading)
		f_error <- y_response - y_est
		
		x_matrix <- e_error
		y_response <- f_error
	}
	
	if (wf == TRUE) {
		###### ADD
		#PLS_core_write_file(w_dist, file = "PLS_w_dist.txt")
		######
		PLS_core_write_file(q_mat, file = "PLS_q_mat.txt")
		PLS_core_write_file(p_mat, file = "PLS_p_mat.txt")
		PLS_core_write_file(w_mat, file = "PLS_w_mat.txt")
		PLS_core_write_file(t_mat, file = paste("PLS_t_mat_learn",file_name,".txt",sep=""))		
	}
	
	return(features_selected)
}

PLS_core_k_learning <- function(x_matrix, y_response, total_f, th, fs=FALSE, wf=TRUE)
{
	n_iter <- 10
	
#	print("evaluating p values")
	pval_vec <- PLS_core_learning(x_matrix,y_response,n_iter,w_adjustment=PLS_adj_top_k_features,k_features=10,get_cor=TRUE)
#	print("evaluating w dist")
	
	w_dist <- c()
	for (i in 1:n_iter) {
		if (pval_vec[i] <= th) {
			w_dist <- c(w_dist,-log(pval_vec[i]))
		}
	}
	
	if (length(w_dist) == 0) {
		print("No component has reached required p-val threshold")
		return(FALSE)
	}
	
	if (length(w_dist) == 1) {
		print("Only one component was found")
	}
	
	
	w_dist <- w_dist / sum(w_dist)
	w_dist <- w_dist * total_f
	w_dist <- round(w_dist)
	w_dist <- w_dist[w_dist > 0]
	
#	print(w_dist)
#	print("calculating w and t")
	res <- PLS_core_k_dist_learning(x_matrix, y_response, w_dist)			
	return(w_dist)
} 

PLS_core_predict <- function(x_vec, num_iterations, wf = FALSE)
{
	p_mat <- PLS_core_read_file("PLS_p_mat.txt")
	w_mat <- PLS_core_read_file("PLS_w_mat.txt")
	t_mat <- rep(0, num_iterations)
	x_mat <- PLS_core_scale_vector(x_vec) #not correct
	
#	print("starting loop")
	
	for (i in 1:num_iterations) {
		#print("t_factor")
		t_factor <- t(x_mat) %*% w_mat[,i]
		t_mat[i] <- t_factor
		#print("new x")
		x_mat <- x_mat - (t_factor %*% t(p_mat[,i]))
		x_mat <- as.vector(x_mat)
		
	}
	
	if (wf == TRUE) {
		PLS_core_write_file(t_mat, file = "PLS_t_mat_predict.txt")
	}
	
}



PLS_core_scale_vector <- function(vec)
{
	res_vec <- vec - mean(vec)
	return(res_vec)
}

PLS_core_scale_matrix <- function(mat)
{
	tmp_mat <- mat
	for (j in 1:ncol(mat)) {
		tmp_col <- PLS_core_scale_vector(tmp_mat[,j])
		tmp_mat[,j] <- tmp_col	
	}
	return(tmp_mat)
}

PLS_core_get_w_loading <- function(x_mat, y_res, scaled = TRUE)
{
	if (scaled == FALSE) {
		x_mat <- PLS_core_scale_matrix(x_mat)
		y_res <- PLS_core_scale_vector(y_res)
	}
	
	c_val <- as.numeric((t(y_res) %*% x_mat %*% t(x_mat) %*% y_res)^(-0.5))
	loading_res <- c_val * (t(x_mat) %*% y_res)
	return(loading_res)
}

PLS_core_write_file <- function(data, file_name, path=paste(Sys.getenv("temp"),"\\",sep="") )
{
	full_name <- paste(path,file_name,sep="")
	write.table(as.matrix(data),file = full_name,sep="\t",eol="\n",row.names=FALSE,col.names=FALSE)
}

PLS_core_read_file <- function(file_name, path=paste(Sys.getenv("temp"),"\\",sep=""))
{
	full_name <- paste(path,file_name,sep="")
	mat <- as.matrix(read.table(full_name,sep="\t"))
	return(mat)
}

PLS_core_default_w_adjustment <- function(w_loading, ...)
{
	return(w_loading)
}

PLS_core_update_features_selected <- function(w_loading, features_selected)
{
	features_selected[w_loading != 0] = features_selected[w_loading != 0] + 1
	return(features_selected)
}

PLS_get_features_by_w <- function(x_mat,learn,num_c,...)
{
	w_mat <- PLS_core_read_file("PLS_w_mat.txt")
#	print("USING top features")
#	print(dim(w_mat))
	f_vector <- rep(FALSE, nrow(w_mat))
#	print(num_c)
	for (i in 1:num_c) {
		tmp_vec <- w_mat[,i]
		tmp_vec[tmp_vec==0] <- FALSE
		tmp_vec[tmp_vec!=0] <- TRUE
		f_vector <- f_vector | tmp_vec
	}
	if (learn == TRUE) {
		new_x_mat <- x_mat[,f_vector]
	}
	else {
		new_x_mat <- x_mat[f_vector]
	}
#	print(dim(new_x_mat))
	return(new_x_mat)
}

PLS_get_t_factors <- function(x_mat,learn,num_c,...)
{
	#### ADD
	#w_dist <- PLS_core_read_file("PLS_w_dist.txt")
	####
	if (learn == TRUE) {
		t_mat <- PLS_core_read_file("PLS_t_mat_learn.txt")
		#### ADD
		#for (i in 1:num_c) {
		#	t_mat[,i] <- t_mat[,i] * (w_dist[i] / sum(w_dist))
		#}
		####
		return(t_mat[,1:num_c])
	}
	else {
		t_mat <- PLS_core_read_file("PLS_t_mat_predict.txt")
		#### ADD
		#for (i in 1:num_c) {
		#	t_mat[i] <- t_mat[i] * (w_dist[i] / sum(w_dist))
		#}
		####	
		return(t_mat[1:num_c])
	}
}



PLS_core_get_loadings <- function(mat, t_factor)
{
	cur_loading <-  PLS_core_calc_loadings_using_LS(mat, t_factor)
	return(cur_loading)
}

PLS_core_calc_loadings_using_LS <- function(org_data, t_factor)
{
	loading_res <- (t(org_data) %*% t_factor) / as.numeric((t(t_factor) %*% t_factor))
	return(loading_res)
}

PLS_find_last_index_vector <- function(num_comp)
{
	w_mat <- PLS_core_read_file("PLS_w_mat.txt")
	len <- nrow(w_mat)
	f_vector <- rep(FALSE, len)
	for (i in 1:num_comp) {
		tmp_vec <- w_mat[,i]
		tmp_vec[tmp_vec==0] <- FALSE
		tmp_vec[tmp_vec!=0] <- TRUE
		f_vector <- f_vector | tmp_vec
	}
	res_vec <- rep(0,len)
	res_vec[f_vector] <- 1
	return(res_vec)
}

####################################################################################################################################
#
# finding top k features using simple sort
#
####################################################################################################################################
PLS_adj_top_k_features <- function(w_loading, k_features, ...)
{
	#HILL
	#res <- PLS_adj_hill_climbing(w_loading,k_features=k_features,...)
	#res <- PLS_adj_sa(w_loading,k_features=k_features,...)
	
	#return(res)
	
	index_vec <- PLS_get_index_vector(w_loading^2, k_features, decreasing=TRUE)
	w_loading[!index_vec] <- 0
	w_loading <- PLS_adj_normalize_vector(w_loading)
	return (w_loading)	
}


PLS_adj_normalize_vector <- function(w_vec)
{
	w_res <- w_vec / sqrt(sum(w_vec^2))
	return(w_res)
}


##
# SVM support functions
##

run_loocv_multi <- function(mat, alg, fvector, ...)
{
	total_error <- 0
	for (i in 1:length(fvector)) {
		#print(i)
		#print("Learning...")
		if (length(dim(mat)) > 0) {
			alg_res <- alg(mat[-i,],fvector[-i],...)
		}
		else {
			alg_res <- alg(mat[-i],fvector[-i],...)
		}
		#print("Predicting...")
		if (length(dim(mat)) > 0) {
			predict_res <- predict(alg_res, t(as.matrix(mat[i,])))
		}
		else {
			predict_res <- predict(alg_res, t(as.matrix(mat[i])))
		}
		if (predict_res != fvector[i]) {
			total_error <- total_error + 1
			#full_text <- paste(as.character(predict_res)," instead of ",as.character(fvector[i]),sep="")
			#print(full_text)
		}
	}
	return(total_error)
}

find_c_regularization <- function(mat, fvector, ...)
{
	#c_grid <- c(0.1,1,10,100,1000,10000)
	c_grid <- c(0.1,1,10)
	
	res_vec <- rep(0,length(c_grid))
	max_index <- 1
	
#	print(dim(mat))
	
	for (i in 1:length(c_grid)) {
#		print(c_grid[i])
		res_vec[i] <- run_loocv_multi(mat, svm, fvector, cost=c_grid[i], ...)
		if (res_vec[i] < res_vec[max_index]) max_index <- i
	}
	
	return(c_grid[max_index])
}


### KNN helper functions

choose_k_for_knn <- function(x_mat, y_fact, k_vec = c(1,3,5,7))
{
	best_k <- 0
	best_score <- 0
	
	for (k in k_vec) {
		res <- knn.cv(x_mat, y_fact, k)
		checks <- res == y_fact
		cur_score <- length(checks[checks])
		if (cur_score > best_score) {
			best_score <- cur_score
			best_k <- k
		} 	
	}
	return(best_k)
}
