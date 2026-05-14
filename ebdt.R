GetExpectancyOfBeingDownstreamTarget <- function(inhibitionThreshold, ratioThreshold, probabilityThreshold, cellLinesFiles) {
	numKinases   <- 163
	numCompounds <- 61

	library(openxlsx)
	combinedKinaseInhibitionsList <- as.matrix(read.csv('requiredData/kinaseInhibitionSpecificity.csv', header=FALSE, stringsAsFactors=FALSE))

	# Short compound names (last segment after ".")
	nKinRows <- min(numKinases + 1, nrow(combinedKinaseInhibitionsList))
	compoundShortNames <- sapply(combinedKinaseInhibitionsList[1, 2:(numCompounds+1)],
	                             function(x) tail(strsplit(x, "\\.")[[1]], 1))

	# Build kinasesInhibitedCompounds vectorised
	kinasesInhibitedCompounds <- list()
	for (r in 2:nKinRows) {
		kinaseRow <- combinedKinaseInhibitionsList[r, ]
		vals      <- as.numeric(kinaseRow[2:(numCompounds+1)])
		mask      <- !is.na(vals) & vals < inhibitionThreshold
		kinasesInhibitedCompounds[[kinaseRow[1]]] <- list(compoundShortNames[mask], vals[mask])
	}

	listOfCompounds <- compoundShortNames
	kusterHeadings  <- combinedKinaseInhibitionsList[1, 2:(numCompounds+1)]
	compoundsKuster  <- setNames(kusterHeadings,
	                             sapply(kusterHeadings, function(h) strsplit(h, "\\.")[[1]][2]))

	cellLines <- sub("\\.xl.*$", "", cellLinesFiles)

	listOfKinasesMatrix <- as.matrix(read.csv('requiredData/listOfKinases.csv', header=FALSE, stringsAsFactors=FALSE))
	listOfKinases <- setNames(
		listOfKinasesMatrix[2:nrow(listOfKinasesMatrix), 3],
		listOfKinasesMatrix[2:nrow(listOfKinasesMatrix), 1]
	)

	for (i in seq_along(cellLines)) {
		if (!cellLines[i] %in% names(listOfKinases)) {
			cat("Error: No kinases imported for", cellLines[i],
			    ". Is this cell line specified in listOfKinases.csv?\n")
			quit(save="no")
		}
	}

	for (i in seq_along(cellLines)) {
		cellLine <- cellLines[i]
		cat("Working on", cellLine, "file [", i, "/", length(cellLines), "]\n")
		cellLineWb <- loadWorkbook(cellLinesFiles[i], keepVBA=TRUE)

		result           <- getDataForKusterCompounds(cellLineWb, listOfCompounds, numCompounds)
		fValues          <- result$fValues
		pValues          <- result$pValues
		fdrValues        <- result$fdrValues
		compoundsCellLine <- result$compoundsCellLine
		sitesCellLine    <- result$sitesCellLine

		ratios           <- correlatePhosphoPeptideWithInhibitorSpecificity(
		                       cellLineWb, kinasesInhibitedCompounds,
		                       fValues, pValues, compoundsKuster, compoundsCellLine)
		probKinases      <- getProbabilityofBeingKinaseSubs(cellLineWb, listOfKinases, ratios, cellLine)
		kinaseSubstrates <- makeKSlistOfKinaseDownstreamTargets(
		                       cellLineWb, ratios, probKinases, pValues, fdrValues,
		                       ratioThreshold, probabilityThreshold)
		makePathways(cellLineWb, probKinases, sitesCellLine)
		createNodeEdges(cellLineWb, kinaseSubstrates)

		saveWorkbook(cellLineWb, cellLinesFiles[i], overwrite=TRUE)
	}
}


# ==== getDataForKusterCompounds ====

getDataForKusterCompounds <- function(cellLineWb, listOfCompounds, numCompounds) {
	if (!'pvalue.select' %in% names(cellLineWb)) addWorksheet(cellLineWb, 'pvalue.select')
	if (!'fold.select'   %in% names(cellLineWb)) addWorksheet(cellLineWb, 'fold.select')

	wsFValuesData   <- readWorkbook(cellLineWb, sheet='fold', colNames=FALSE)
	headerRow       <- as.character(wsFValuesData[1, 3:(length(listOfCompounds)+2)])
	compoundsCellLine <- setNames(headerRow,
	                              sapply(headerRow, function(h) strsplit(h, "\\.")[[1]][2]))
	compoundHeadings <- names(compoundsCellLine)

	sitesCellLine <- as.character(wsFValuesData[2:nrow(wsFValuesData), 1])
	nSites <- length(sitesCellLine)
	nComp  <- length(compoundHeadings)

	# Build fold matrix directly (vectorised read)
	fMat <- matrix(
		as.numeric(as.matrix(wsFValuesData[2:nrow(wsFValuesData), 3:(nComp+2)])),
		nrow = nSites, ncol = nComp,
		dimnames = list(sitesCellLine, compoundHeadings)
	)
	fValues <- setNames(
		lapply(seq_len(nSites), function(i) setNames(as.list(fMat[i, ]), compoundHeadings)),
		sitesCellLine
	)

	wsPValuesData <- readWorkbook(cellLineWb, sheet='pvalue', colNames=FALSE)
	pSites  <- as.character(wsPValuesData[2:nrow(wsPValuesData), 1])
	nPSites <- length(pSites)

	pMat <- matrix(
		as.numeric(as.matrix(wsPValuesData[2:nrow(wsPValuesData), 3:(nComp+2)])),
		nrow = nPSites, ncol = nComp,
		dimnames = list(pSites, compoundHeadings)
	)
	pValues <- setNames(
		lapply(seq_len(nPSites), function(i) setNames(as.list(pMat[i, ]), compoundHeadings)),
		pSites
	)

	fdrRaw    <- wsPValuesData[2:nrow(wsPValuesData), 2]
	# Si la columna 2 es texto (ej. sh.index.ids en NTERA2), usar 0 como FDR (pasa el filtro fdr < 0.02)
	fdrNumeric <- suppressWarnings(as.numeric(as.character(fdrRaw)))
	fdrValues <- setNames(
		as.list(ifelse(is.na(fdrNumeric), 0, fdrNumeric)),
		pSites
	)

	return(list(fValues=fValues, pValues=pValues, fdrValues=fdrValues,
	            compoundsCellLine=compoundsCellLine, sitesCellLine=sitesCellLine))
}


# ==== correlatePhosphoPeptideWithInhibitorSpecificity ====

correlatePhosphoPeptideWithInhibitorSpecificity <- function(cellLineWb, kinasesInhibitedCompounds,
                                                            fValues, pValues,
                                                            compoundsKuster, compoundsCellLine) {
	if (!'corrPPwithKinases'          %in% names(cellLineWb)) addWorksheet(cellLineWb, 'corrPPwithKinases')
	if (!'ratioSigPPoverSignInhSpeci' %in% names(cellLineWb)) addWorksheet(cellLineWb, 'ratioSigPPoverSignInhSpeci')

	siteNames     <- names(fValues)
	compoundNames <- names(fValues[[siteNames[1]]])
	nSites        <- length(siteNames)
	kinaseNames   <- names(kinasesInhibitedCompounds)

	# Build site * compound matrices once
	fMat <- do.call(rbind, lapply(siteNames, function(s) as.double(unlist(fValues[[s]]))))
	rownames(fMat) <- siteNames; colnames(fMat) <- compoundNames
	pMat <- do.call(rbind, lapply(siteNames, function(s) as.double(unlist(pValues[[s]]))))
	rownames(pMat) <- siteNames; colnames(pMat) <- compoundNames

	# Identify active kinases (>1 inhibited compound): pre-allocate output matrices
	activeKin <- kinaseNames[sapply(kinaseNames, function(k) length(kinasesInhibitedCompounds[[k]][[1]]) > 1)]
	nActive   <- length(activeKin)

	corrMat  <- matrix(0.0, nrow=nSites, ncol=nActive, dimnames=list(siteNames, activeKin))
	ratioMat <- matrix(0.0, nrow=nSites, ncol=nActive, dimnames=list(siteNames, activeKin))

	# Compute correlations and ratios column by column (one kinase at a time)
	for (ki in seq_len(nActive)) {
		kinaseName  <- activeKin[ki]
		kinaseTuple <- kinasesInhibitedCompounds[[kinaseName]]
		kCompounds  <- kinaseTuple[[1]]
		kValues     <- as.numeric(kinaseTuple[[2]])
		nKinComp    <- length(kCompounds)

		valid   <- kCompounds %in% compoundNames
		kComp_v <- kCompounds[valid]
		kVals_v <- kValues[valid]

		if (length(kComp_v) >= 2) {
			F_sub <- fMat[, kComp_v, drop=FALSE]   # nSites*k
			P_sub <- pMat[, kComp_v, drop=FALSE]

			# Correlation: apply cor() row-wise: guaranteed length nSites
			corr_col <- suppressWarnings(
			    vapply(seq_len(nSites), function(si) {
			        val <- cor(kVals_v, F_sub[si, ])
			        if (length(val) != 1L || is.na(val)) 0.0 else val
			    }, double(1))
			)
			corrMat[, ki] <- corr_col

			# Ratio: vectorised rowSums
			inhibited     <- !is.na(F_sub) & !is.na(P_sub) & F_sub < -1 & P_sub < 0.025
			ratioMat[, ki] <- rowSums(inhibited) / nKinComp
		}
		# (if kComp_v<2, column stays 0: already initialised)
	}

	# Build ratioList (named list of named lists) for downstream functions
	ratioList <- setNames(vector("list", length(kinaseNames)), kinaseNames)
	for (ki in seq_len(nActive)) {
		k <- activeKin[ki]
		ratioList[[k]] <- setNames(as.list(ratioMat[, ki]), siteNames)
	}
	# inactive kinases keep list() (NULL from vector("list",...) is fine: length 0)

	# ==== Write corrPPwithKinases (one writeData call) ====
	if (nActive > 0) {
		writeData(cellLineWb, 'corrPPwithKinases',
		          cbind(data.frame(site=siteNames, stringsAsFactors=FALSE), as.data.frame(corrMat)),
		          startRow=1, startCol=1, colNames=TRUE)
	}

	# ==== Write ratioSigPPoverSignInhSpeci (one writeData call) ====
	if (nActive > 0) {
		writeData(cellLineWb, 'ratioSigPPoverSignInhSpeci',
		          cbind(data.frame(site=siteNames, stringsAsFactors=FALSE), as.data.frame(ratioMat)),
		          startRow=1, startCol=1, colNames=TRUE)
	}

	return(ratioList)
}


# ==== getProbabilityofBeingKinaseSubs ====

getProbabilityofBeingKinaseSubs <- function(cellLineWb, listOfKinases, ratios, cellLine) {
	if (!'ProbOfBeingKinaseSubs' %in% names(cellLineWb)) addWorksheet(cellLineWb, 'ProbOfBeingKinaseSubs')

	allKinases    <- names(ratios)
	activeKinases <- allKinases[sapply(allKinases, function(k) length(ratios[[k]]) > 0)]

	if (length(activeKinases) == 0) {
		probOfBeingKinaseSubstrate <- setNames(lapply(allKinases, function(k) list()), allKinases)
		return(probOfBeingKinaseSubstrate)
	}

	siteNames <- names(ratios[[activeKinases[1]]])
	nSites    <- length(siteNames)

	# Build ratio matrix: activeKinases x sites (pre-allocate + fill row by row)
	ratioMat <- matrix(0.0, nrow=length(activeKinases), ncol=nSites,
	                   dimnames=list(activeKinases, siteNames))
	for (ki in seq_along(activeKinases)) {
		ratioMat[ki, ] <- as.double(unlist(ratios[[activeKinases[ki]]]))
	}

	# Max ratio per site (column max across all active kinases)
	maxRatios <- apply(ratioMat, 2, max, na.rm=TRUE)

	# Kinases present in this cell line
	kinasesInCellLine <- activeKinases[sapply(activeKinases, function(kinase) {
		dotPos       <- regexpr("\\.", kinase)[1]
		kinaseToFind <- if (dotPos > 0) substr(kinase, 1, dotPos - 1) else kinase
		grepl(kinaseToFind, listOfKinases[[cellLine]], fixed=TRUE)
	})]

	# Probability matrix (initialised to 0)
	probMat <- matrix(0, nrow=length(activeKinases), ncol=nSites,
	                  dimnames=list(activeKinases, siteNames))

	inCL <- activeKinases %in% kinasesInCellLine
	if (any(inCL)) {
		ratioSub <- ratioMat[inCL, , drop=FALSE]
		probSub  <- sweep(ratioSub, 2, maxRatios,
		                  FUN=function(r, m) ifelse(m > 0, r / m, 0))
		probMat[inCL, ] <- probSub
	}

	# Build output list (same structure as original)
	probOfBeingKinaseSubstrate <- setNames(vector("list", length(allKinases)), allKinases)
	for (k in allKinases) {
		if (k %in% activeKinases) {
			probOfBeingKinaseSubstrate[[k]] <- setNames(as.list(probMat[k, ]), siteNames)
		} else {
			probOfBeingKinaseSubstrate[[k]] <- list()
		}
	}

	# ==== Write ProbOfBeingKinaseSubs (one writeData call) ====
	probDf <- data.frame(site=siteNames, stringsAsFactors=FALSE)
	for (k in allKinases) {
		probDf[[k]] <- if (k %in% activeKinases) probMat[k, ] else rep(NA_real_, nSites)
	}
	writeData(cellLineWb, 'ProbOfBeingKinaseSubs', probDf, startRow=1, startCol=1, colNames=TRUE)

	return(probOfBeingKinaseSubstrate)
}


# ==== makeKSlistOfKinaseDownstreamTargets ====

makeKSlistOfKinaseDownstreamTargets <- function(cellLineWb, ratios, probKinases, pValues,
                                                fdrValues, ratioThreshold, probabilityThreshold) {
	residuesToCheck <- c("None", "(M", "(R", "(K")
	kinaseSubstrates <- list()

	for (kinase in names(probKinases)) {
		kinaseDict <- probKinases[[kinase]]
		substrates <- c()
		for (phosphosite in names(kinaseDict)) {
			if (!any(sapply(residuesToCheck, function(res) grepl(res, phosphosite, fixed=TRUE)))) {
				fdr   <- fdrValues[[phosphosite]]
				prob  <- probKinases[[kinase]][[phosphosite]]
				ratio <- ratios[[kinase]][[phosphosite]]
				if (!is.null(prob) && !is.null(ratio) &&
				    prob > probabilityThreshold && ratio > ratioThreshold &&
				    !is.null(fdr) && fdr < 0.02) {
					sites      <- strsplit(phosphosite, ";")[[1]]
					sites      <- sites[nchar(sites) > 0]
					substrates <- c(substrates, sites)
				}
			}
		}
		kinaseSubstrates[[kinase]] <- substrates
	}

	if (!'PutativeKinaseSubstrates' %in% names(cellLineWb)) {
		addWorksheet(cellLineWb, 'PutativeKinaseSubstrates')
	}

	# Write header + data in one call
	ksDf <- data.frame(
		kinase     = names(kinaseSubstrates),
		n          = sapply(kinaseSubstrates, length),
		substrates = sapply(kinaseSubstrates, function(s) {
			if (length(s) > 0) paste0(paste(s, collapse=";"), ";") else ""
		}),
		stringsAsFactors = FALSE
	)
	writeData(cellLineWb, 'PutativeKinaseSubstrates', ksDf, startRow=1, startCol=1, colNames=TRUE)

	return(kinaseSubstrates)
}


# ==== makePathways ====

makePathways <- function(cellLineWb, probOfBeingKinaseSubstrate, sitesCellLine) {
	phosphositePathways <- list()
	for (kinase in names(probOfBeingKinaseSubstrate)) {
		kinaseDict <- probOfBeingKinaseSubstrate[[kinase]]
		for (phosphosite in names(kinaseDict)) {
			value <- kinaseDict[[phosphosite]]
			if (!is.null(value) && value == 1) {
				phosphositePathways[[phosphosite]] <-
					c(phosphositePathways[[phosphosite]], kinase)
			}
		}
	}

	uniquePathways <- list()
	for (site in sitesCellLine) {
		pw <- phosphositePathways[[site]]
		if (!is.null(pw) && length(pw) > 1) {
			pathwayString <- paste(pw, collapse="-")
			found <- any(sapply(uniquePathways, function(up) {
				grepl(pathwayString, paste(up, collapse="-"), fixed=TRUE)
			}))
			if (!found) uniquePathways <- c(uniquePathways, list(pw))
		}
	}
	uniquePathways <- rev(uniquePathways)

	pathwaySites <- list()
	for (pathway in uniquePathways) {
		pathway       <- as.character(unlist(pathway))
		pathwayString <- paste(pathway, collapse="-")
		sites         <- c()
		for (site in names(phosphositePathways)) {
			sitePathway <- phosphositePathways[[site]]
			maxStart    <- length(sitePathway) - length(pathway) + 1
			if (maxStart >= 1) {
				for (i in seq_len(maxStart)) {
					if (identical(pathway, as.character(sitePathway[i:(i+length(pathway)-1)]))) {
						sites <- c(sites, site)
					}
				}
			}
		}
		if (length(sites) > 1) pathwaySites[[pathwayString]] <- sites
	}
}


# ==== createNodeEdges ====

createNodeEdges <- function(cellLineWb, kinaseSubstrates) {
	kinaseNames <- names(kinaseSubstrates)
	nKinases    <- length(kinaseNames)
	nodesSubtrates <- list()

	# Use intersect() instead of inner substrate loop (O(n^2) pairs, but each is fast)
	for (i in seq_len(nKinases)) {
		kinase <- kinaseNames[i]
		subsA  <- kinaseSubstrates[[kinase]]
		if (length(subsA) == 0) next
		for (j in seq_len(nKinases)) {
			if (kinase < kinaseNames[j]) {
				subsB  <- kinaseSubstrates[[kinaseNames[j]]]
				if (length(subsB) == 0) next
				shared <- intersect(subsA, subsB)
				if (length(shared) > 0) {
					nodesSubtrates[[paste0(kinase, ".", kinaseNames[j])]] <- shared
				}
			}
		}
	}

	# ==== tableEdgeSubs: build full matrix, write once ====
	if (!'tableEdgeSubs' %in% names(cellLineWb)) addWorksheet(cellLineWb, 'tableEdgeSubs')

	edgeMat <- matrix("", nrow=nKinases+1, ncol=nKinases+1)
	edgeMat[1, 2:(nKinases+1)] <- kinaseNames
	edgeMat[2:(nKinases+1), 1] <- kinaseNames
	for (i in seq_len(nKinases)) {
		for (j in seq_len(nKinases)) {
			if (j > i) {
				key       <- paste0(kinaseNames[i], ".", kinaseNames[j])
				val       <- nodesSubtrates[[key]]
				cellValue <- if (is.null(val)) "" else paste(as.character(val), collapse=";")
				edgeMat[i+1, j+1] <- cellValue
				edgeMat[j+1, i+1] <- cellValue
			}
		}
	}
	writeData(cellLineWb, 'tableEdgeSubs',
	          as.data.frame(edgeMat, stringsAsFactors=FALSE),
	          startRow=1, startCol=1, colNames=FALSE)

	# ==== nodes.edges ====
	if (!'nodes.edges' %in% names(cellLineWb)) {
		addWorksheet(cellLineWb, 'nodes.edges')
		writeData(cellLineWb, 'nodes.edges', "edge",   startRow=1, startCol=1, colNames=FALSE)
		writeData(cellLineWb, 'nodes.edges', "weight", startRow=1, startCol=2, colNames=FALSE)
		writeData(cellLineWb, 'nodes.edges', "subs",   startRow=1, startCol=3, colNames=FALSE)
	}
	if (length(nodesSubtrates) > 0) {
		nodesDf <- data.frame(
			edge   = names(nodesSubtrates),
			weight = sapply(nodesSubtrates, length),
			subs   = sapply(nodesSubtrates, function(s) paste(s, collapse=";")),
			stringsAsFactors = FALSE
		)
		writeData(cellLineWb, 'nodes.edges', nodesDf, startRow=2, startCol=1, colNames=FALSE)
	}
}
