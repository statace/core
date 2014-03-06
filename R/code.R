sa.corr <- function(..., pearson = TRUE, kendall = TRUE, spearman = TRUE, one.tailed = FALSE) {
	vars <- list(...)
	names(vars) <- sa.get.dot.parameter.names(...)

	# Ensure the variables are of the same length
	lengths <- sapply(vars, length)
	if (min(lengths) < 1 || min(lengths) != max(lengths)) {
		cat("The variables need to be of the same length")
		return(NULL)
	}

	x <- as.data.frame(vars)
	num.vars <- length(vars)

	pearson = as.logical(pearson)
	kendall = as.logical(kendall)
	spearman = as.logical(spearman)
	one.tailed = as.logical(one.tailed)

	if (pearson[1] == FALSE && kendall[1] == FALSE && spearman[1] == FALSE) {
		return(NULL)
	}

	# The number of input variables
	num.vars = length(x)

	# Output is going to be a matrix
	output = numeric(0)

	if (pearson[1] == TRUE) {
		result = corr.test(x = x, method = "pearson", adjust = "none")
		r = sa.add.matrix.row.subname(sa.format.all.to.all.matrix(result$r), "Pearson's r")
		p = sa.add.matrix.row.subname(sa.format.all.to.all.matrix(result$p), "Sig. of Pearson's r")
		if (one.tailed[1] == TRUE) {
			p = p / 2
		}
		output = r
		output = sa.vertical.matrix.merge(output, p)
	}

	if (kendall[1] == TRUE) {
		result = corr.test(x = x, method = "kendall", adjust = "none")
		r = sa.add.matrix.row.subname(sa.format.all.to.all.matrix(result$r), "Kendall's tau")
		p = sa.add.matrix.row.subname(sa.format.all.to.all.matrix(result$p), "Sig. of Kendall's tau")
		if (one.tailed[1] == TRUE) {
			p = p / 2
		}
		if (length(output) > 0) {
			output = sa.vertical.matrix.merge(output, r)
		}
		else {
			output = r
		}
		output = sa.vertical.matrix.merge(output, p)
	}

	if (spearman[1] == TRUE) {
		result = corr.test(x = x, method = "spearman", adjust = "none")
		r = sa.add.matrix.row.subname(sa.format.all.to.all.matrix(result$r), "Spearman's rho")
		p = sa.add.matrix.row.subname(sa.format.all.to.all.matrix(result$p), "Sig. of Spearman's rho")
		if (one.tailed[1] == TRUE) {
			p = p / 2
		}
		if (length(output) > 0) {
			output = sa.vertical.matrix.merge(output, r)
		}
		else {
			output = r
		}
		output = sa.vertical.matrix.merge(output, p)
	}

	# Generate a counts matrix from the single count number (cols & rows are the # of variables)
	n = matrix(result$n, nrow = num.vars, ncol = num.vars, dimnames = list(names(x)))
	n = sa.format.all.to.all.matrix(n)
	n = sa.add.matrix.row.subname(n, "Count")

	# Merge the counts matrix into the result
	output = sa.vertical.matrix.merge(output, n)

	# Convert matrix to data frame, and add class
	output <- as.data.frame(output)

	# Mark as a StatAce object
	output <- sa.mark.object(output)

	# Set data frame to not highlight rows with NAs
	output <- sa.set.no.na.row.highlight(output)

	# Set data frame to not print NAs
	output <- sa.set.no.na.print(output)

	output
}

sa.format.all.to.all.matrix <- function(x) {
	x = as.matrix(x)
	diag(x) <- NA
	x
}

sa.desc <- function(..., freq.table = F, chart = "") {
	vars <- list(...)
	var.names <- sa.get.dot.parameter.names(...)
	num.vars <- length(vars)

	count            <- rep(NA_integer_, num.vars)
	valid            <- rep(NA_integer_, num.vars)
	missing          <- rep(NA_integer_, num.vars)
	mean             <- rep(NA_real_, num.vars)
	sd               <- rep(NA_real_, num.vars)
	variance         <- rep(NA_real_, num.vars)
	med              <- rep(NA_real_, num.vars)
	mode             <- rep(NA_character_, num.vars)
	min              <- rep(NA_real_, num.vars)
	max              <- rep(NA_real_, num.vars)
	range            <- rep(NA_real_, num.vars)
	quartile.1       <- rep(NA_real_, num.vars)
	quartile.3       <- rep(NA_real_, num.vars)
	std.err.of.mean  <- rep(NA_real_, num.vars)

	if (! chart %in% c("pie", "bar", "hist", "hist.norm.curve")) {
		chart = NULL
	}
	
	i <- 1

	for (i in 1:num.vars) {
		# Get the variable
		var <- vars[[i]]

		# Get counts
		count[i] <- length(var)

		if (is.vector(var)) {
			missing[i] <- sum(is.na(var))
			valid[i] <- count[i] - missing[i]

			# Determine the mode
			unique <- unique(var)
			tabulated <- tabulate(match(var, unique))
			if (max(tabulated) > 1) {
				mode[i] <- as.character(unique[which.max(tabulated)])
			}
		
			# All stats except the mode are only applicable if the variable is numeric
			if (is.numeric(var)) {
				# Get the summary from the psych package
				psych <- as.list(describe(var))
	
				mean[i] <- psych[["mean"]]
				sd[i] <- psych[["sd"]]
				variance[i] <- sd[i] ^ 2
				med[i] <- psych[["median"]]
				min[i] <- psych[["min"]]
				max[i] <- psych[["max"]]
				range[i] <- psych[["range"]]

				base <- as.list(summary(var))
				quartile.1[i] <- base[["1st Qu."]]
				quartile.3[i] <- base[["3rd Qu."]]

				std.err.of.mean[i] <- psych[["se"]]
			}
		}
	}

	# Create descriptives data frame
	desc.df <- data.frame(Count_Total = count, Count_Valid = valid, Count_Missing = missing, Mean = mean,
		Std_Dev = sd, Variance = variance, Median = med, Mode = mode, Min = min, Max = max, Range = range,
		Quartile_1 = quartile.1, Quartile_3 = quartile.3, Std_Err_of_Mean = std.err.of.mean,
		row.names = var.names)

	# Set data frame as transposed
	desc.df <- sa.set.transposed(desc.df)

	# Set data frame to not highlight rows with NAs
	desc.df <- sa.set.no.na.row.highlight(desc.df)

	# Set data frame to not print NAs
	desc.df <- sa.set.no.na.print(desc.df)

	# If no frequency table and no chart, return the descriptives data frame directly
	if (!freq.table && is.null(chart)) {
		# Mark as a StatAce object
		desc.df <- sa.mark.object(desc.df)

		desc.df
	}
	else {
		# There should be a frequency table and/or a chart -> prepare a list to return
		result <- list()
		result[["Descriptives"]] <- desc.df

		if (freq.table) {
			# Add a frequency table for every variable to the result list
			for (i in 1:num.vars) {
				result[[paste("Frequency table: ", var.names[i], sep = "")]] <- sa.freq.table(vars[[i]])
			}
		}

		if (!is.null(chart)) {
			# There should be a chart

			for (i in 1:num.vars) {
				var <- vars[[i]]
				name <- var.names[i]

				if (is.vector(var)) {
					if (chart == "pie") {
						# Print the chart
						sa.freq.pie(var, name);

						# Add a placeholder for the plot
						result[[paste("Pie chart:", name)]] <- sa.place.plot()
					}
					# In case the var is not numeric, always stop with a bar chart
					else if (chart == "bar" || !is.numeric(var)) {
						# Print the chart
						sa.freq.bar(var, name);

						# Add a placeholder for the plot
						result[[paste("Bar chart:", name)]] <- sa.place.plot()
					}
					else if (chart == "hist") {
						# Print the chart
						sa.freq.hist(var, name);

						# Add a placeholder for the plot
						result[[paste("Histogram:", name)]] <- sa.place.plot()
					}
					else if (chart == "hist.norm.curve") {
						# Print the chart
						sa.freq.hist.norm.curve(var, name);

						# Add a placeholder for the plot
						result[[paste("Histogram with normal curve:", name)]] <- sa.place.plot()
					}
				}
			}
		}

		# Mark as a StatAce object
		result <- sa.mark.object(result)

		result
	}
}

sa.freq.pie <- function(var, title) {
	slices <- table(var)
	names <- names(slices)
	pct <- round(slices / sum(slices) * 100)
	labels <- paste(names, pct)
	labels <- paste(labels, "%", sep = "")
	
	par(mfrow = c(1, 1))
	pie(slices, labels = labels, col = rainbow(length(labels)), main=title)
}

sa.freq.bar <- function(var, title) {
	plot <- ggplot(as.data.frame(var), aes(sa.as.formatted.factor(var))) + geom_bar() + xlab("Value") + ylab("Frequency") + ggtitle(title) + scale_y_continuous(expand = c(0, 0))
	print(plot)
}

sa.freq.hist <- function(var, title) {
	plot <- ggplot(data.frame(var), aes(var)) + geom_histogram(binwidth = (max(var, na.rm = TRUE) - min(var, na.rm = TRUE)) / 20) + xlab("Value") + ylab("Frequency") + ggtitle(title) + scale_y_continuous(expand = c(0, 0))
	print(plot)
}

sa.freq.hist.norm.curve <- function(var, title) {
	binwidth <- (max(var, na.rm = TRUE) - min(var, na.rm = TRUE)) / 20
	plot <- ggplot(data.frame(var), aes(var)) + geom_histogram(binwidth = binwidth) + stat_function(fun = sa.dnorm, args = list(mean = mean(var, na.rm = TRUE), sd = sd(var, na.rm = TRUE), multiplier = length(var) * binwidth)) + xlab("Value") + ylab("Frequency") + ggtitle(title) + scale_y_continuous(expand = c(0, 0))
	print(plot)
}

sa.dnorm <- function (..., mean, sd, multiplier) {
	dnorm(..., mean = mean, sd = sd) * multiplier
}

sa.freq.table <- function(var) {
	missing <- sum(is.na(var))

	# Create a dataframe with the frequency of the values. The Value column is a factor.
	freq.df <- as.data.frame(table(var, dnn = "Value"))

	# Rename the frequency heading to "Frequency"
	names(freq.df)[names(freq.df) == "Freq"] <- "Frequency"

	# Add a column with the % of the valid values
	freq.df <- transform(freq.df, Of.Valid = prop.table(freq.df$Frequency))

	# Add a column with the cumulative % of the valid values
	freq.df <- transform(freq.df, Cumulative.Of.Valid = cumsum(freq.df$Of.Valid))

	# If there are missing values, add a row for them
	if(missing > 0) {
		# Create a data frame for the missing row
		missing.row = data.frame(Value = "[Missing]", Frequency = missing, Of.Valid = NA, Cumulative.Of.Valid = NA)

		# Add the [Missing] string to the factor of the value column in the frequency data frame
		levels(freq.df$Value) <- c(levels(freq.df$Value), "[Missing]")

		# Add the missing row to the end of the frequency data frame
		freq.df <- rbind(freq.df, missing.row)
	}

	# Add a column with the % of all values
	freq.df <- transform(freq.df, Of.All = prop.table(freq.df$Frequency))

	# Reorder the table columns
	freq.df <- freq.df[c("Value", "Frequency", "Of.All", "Of.Valid", "Cumulative.Of.Valid")]

	freq.df
}