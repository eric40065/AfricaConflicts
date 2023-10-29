library(tsne); library(hash); library(kernlab); library(data.table); library(lubridate); library(pracma); library(clue); library(caret); library(maps); library(mapdata)
removeSpeicalCharacters = function(string, Characters = c(" ", "'")){
    for(char in Characters){
        string = gsub(char, "", string)
    }
    return(string)
}
addSoleClassToActors = function(actor){
    actor[actor == ""] = "Solo"
    return(actor)
}
preprocessActor = function(data){
    data$event1_actor1 = addSoleClassToActors(removeSpeicalCharacters(data$event1_actor1))
    data$event1_actor2 = addSoleClassToActors(removeSpeicalCharacters(data$event1_actor2))
    data$event2_actor1 = addSoleClassToActors(removeSpeicalCharacters(data$event2_actor1))
    data$event2_actor2 = addSoleClassToActors(removeSpeicalCharacters(data$event2_actor2))
    # fwrite(data, "/Users/eric/Desktop/UCD/Courses/STA\ 260/Africa_data_7M.csv", row.names = FALSE)
    return(data)
}
createSentence = function(data){
    unique_event1 = unique(data$event1_id)
    row_end = cumsum(table(factor(data$event1_id, levels = unique_event1)))
    row_start = c(1, row_end + 1)
    actor11 = data$event1_actor1
    actor12 = data$event1_actor2
    actor21 = data$event2_actor1
    actor22 = data$event2_actor2
    sentences = lapply(1:length(unique_event1), function(i){
        row_now = row_start[i] : row_end[i]
        paste(paste(actor11[row_now], actor12[row_now], collapse = " "), paste(actor21[row_now], actor22[row_now], collapse = " "), collapse = " ")
    })
    return(sentences)
}
formSentenceFromFile = function(inputPath = NULL, outputPath, data = NULL){
    if(is.null(data) & is.null(inputPath)){error("Missing values. At least one of inputs on inputPath or data is required.")}
    if(is.null(data)){data = as.data.frame(fread(inputPath))}
    ## preprocess the data if there's space in variable actor
    if(sum(grepl(" ", data$event2_actor1[1:min(1000, dim(data)[2])])) > 0){
        pracma::fprintf("Removing space from actor...")
        data = preprocessActor(data)
        pracma::fprintf(" finished. \n")
    }
    pracma::fprintf("Constructing sentence...")
    sentences = createSentence(data)
    write.table(unlist(sentences), file = outputPath, row.names = FALSE, col.names = FALSE, quote = FALSE)
    pracma::fprintf(" finished. \n")
    return(unlist(sentences))
}
CreateNeighborhood = function(inputPath, outputPath = NULL, distance = 30, timeSpan = 30){
    data = as.data.frame(fread(inputPath))
    ## change the event date to yyyy-mm-dd format
    data$event_date = parse_date_time(data$event_date, orders = "d B y")
    
    ## sort the data by event date
    data = data[order(data$event_date), ]
    
    ## create dateDiff: indicate the time difference between a event and the earliest event in the data
    data$dateDiff = as.numeric(difftime(data$event_date, data$event_date[1], units = "days"))
    
    ## rearrange the data
    importantVariables = c("longitude", "latitude", "dateDiff", "actor1", "actor2", "event_date", "data_id")
    if(sum(!(importantVariables %in% colnames(data))) > 0){
        warning("No enough coumns", immediate. = TRUE)
    }
    data = data[, c(importantVariables, colnames(data)[!(colnames(data) %in% importantVariables)])]
    
    ## get the unique start and end days for each event
    date_index_start = unique(approx(data$dateDiff, 1:dim(data)[1], xout = data$dateDiff, ties = "min")$y)
    date_index_end = c(unique(round(approx(data$dateDiff, 1:dim(data)[1], xout = data$dateDiff + timeSpan, ties = "max", rule = 2)$y)), 
                       rep(dim(data)[1], length(unique(data$dateDiff[data$dateDiff + timeSpan > max(data$dateDiff)]))))
    
    ## get pairwise data
    data_location_time = as.matrix(data[, 1:3]) # easy access
    pracma::fprintf("Pairing events...")
    pairwiseData = do.call(rbind.data.frame, lapply(1:length(date_index_start), function(i){
        data_now = data_location_time[date_index_start[i]:date_index_end[i], , drop = FALSE]
        reference = data_now[, 3] == data_now[1, 3]
        distMatrix = rdist::cdist(data_now[, 1:2, drop = FALSE], data_now[reference, 1:2, drop = FALSE]) < distance/111
        # distMatrix = outer(1:sum(data_now[, 3] == data_now[1, 3]), 1:dim(data_now)[1], Vectorize(function(i, j){distm(data_now[i, 1:2], data_now[j, 1:2])/1000})) < distance # exact distance, take plenty of time, don't recommend
        
        repeats = colSums(distMatrix)
        data_now = data[date_index_start[i]:date_index_end[i], , drop = FALSE]
        augmentedData = data.frame(matrix(nrow = sum(distMatrix), ncol = 8))
        event2RowIndex = unlist(apply(distMatrix, 2, which))
        augmentedData$X1 = rep(data_now$data_id[reference], repeats) # event 1 id
        augmentedData$X2 = data_now$data_id[event2RowIndex] # event 2 id
        augmentedData$X3 = data_now$dateDiff[1] # event 1 start day
        augmentedData$X4 = data_now$dateDiff[event2RowIndex] # event 2 start day
        augmentedData$X5 = rep(data_now$actor1[reference], repeats) ## event 1 actor1
        augmentedData$X6 = rep(data_now$actor2[reference], repeats) ## event 1 actor1
        augmentedData$X7 = data_now$actor1[event2RowIndex] # event 2 actor 1
        augmentedData$X8 = data_now$actor2[event2RowIndex] # event 2 actor 2
        return(augmentedData)
    }))
    colnames(pairwiseData) = c("event1_id", "event2_id", "event1_start_date", "event2_start_date", "event1_actor1", "event1_actor2", "event2_actor1", "event2_actor2")
    pracma::fprintf(" finished.\n")
    if(!is.null(outputPath)){
        pracma::fprintf("Writing the file...")
        fwrite(pairwiseData, outputPath, row.names = FALSE)
        pracma::fprintf("finished.\n The file is stored at %s \n", paste(getwd(), outputPath, sep = "/"))
    }
    return(pairwiseData)
}
vec2tsne = function(input_file, output_file, tsne_dim = 2, max_iter = 500, ...){
    embedding = read.table(input_file, quote="\"", comment.char="")
    actorLevels = sort(embedding[, 1], index.return = TRUE)
    embedding = embedding[actorLevels$ix, 2:dim(embedding)[2]]
    actorLevels = actorLevels$x
    tsne_fit = tsne(embedding, k = tsne_dim, max_iter = max_iter, ...)
    return_df = data.frame(actors = actorLevels, tsne_fit)
    colnames(return_df) = c("actors", paste0("Embedding", 1:tsne_dim))
    write.csv(return_df, file = output_file, row.names = FALSE)
    return(return_df)
}
matchVaribles = function(originalDataPath, augmentedDataPath, actorTSNE){
    ##### match actor embedding and create duration
    largeData = as.data.frame(fread(augmentedDataPath))
    data = as.data.frame(fread(originalDataPath))
    data = matchEmbedding(data, actorTSNE)
    data = matchDuration(data, largeData)
    return(data)
}
getCountryData = function(data, country = "Somalia", getEventAfter = NULL, getFatalitiesLargerThan = 0){
    fatalitiesIndicator = (data$fatalities >= getFatalitiesLargerThan)
    if(!all(country %in% data$country)){warning("Some countries specificed are not found in the input data.")}
    countryIndicator = (data$country %in% country)
    countryData = data[countryIndicator & fatalitiesIndicator, c(which(colnames(data) %in% c("fatalities", "latitude", "longitude", "event_date", "actor1", "actor2")), dim(data)[2] - ((2 * dim(actorTSNE)[2] - 2):0))]
    durationIndicator = countryData$duration > 0
    eventDate = parse_date_time(data$event_date[countryIndicator & fatalitiesIndicator], orders = "d B y") ## change the event date to yyyy-mm-dd format
    if(is.null(getEventAfter)){
        getEventAfter = min(eventDate)
    }else{
        getEventAfter = parse_date_time(getEventAfter, orders = "d B y")
    }
    dateIndicator = (eventDate >= getEventAfter)
    countryData = countryData[dateIndicator & durationIndicator, ]
    return(countryData)
}
matchEmbedding = function(data, actorTSNE){
    ## assign the first column in actorTSNE as the names of rows and remove the first column
    row.names(actorTSNE) = removeSpeicalCharacters(tolower(actorTSNE$actors))
    actorTSNE = actorTSNE[, 2:dim(actorTSNE)[2]]
    EmbeddingDim = dim(actorTSNE)[2]
    
    ## select the rows where the actors involved are embedded
    data$actor2 = addSoleClassToActors(data$actor2)
    allActors = unique(c(data$actor1, data$actor2))
    allActors = data.frame(deCase = removeSpeicalCharacters(tolower(allActors)), actors = allActors)
    allActors = allActors[allActors$deCase %in% row.names(actorTSNE), ]
    data = data[data$actor1 %in% allActors$actors & data$actor2 %in% allActors$actors, ]
    
    ## Embed actor1 and actor2
    actor1 = removeSpeicalCharacters(tolower(data$actor1))
    actor2 = removeSpeicalCharacters(tolower(data$actor2))
    data = cbind(data, matrix(0, nrow = dim(data)[1], ncol = EmbeddingDim * 2))
    data[, dim(data)[2] - (2 * EmbeddingDim - 1):EmbeddingDim] = actorTSNE[actor1, ]
    data[, dim(data)[2] - (EmbeddingDim - 1):0] = actorTSNE[actor2, ]
    colnames(data)[dim(data)[2] - (2 * EmbeddingDim - 1):0] = c(paste0("actor1Embedding", 1:EmbeddingDim), paste0("actor2Embedding", 1:EmbeddingDim))
    return(data)
}
matchDuration = function(data, largeData){
    data$duration = sapply(split(largeData$event2_start_date - largeData$event1_start_date, factor(largeData$event1_id, levels = unique(largeData$event1_id))), function(x){median(x[x > 0])})
    data$duration[is.na(data$duration)] = 0
    return(data)
}
myNormalize = function(data, importance = c(2, 1, 1)){
    # data = somalia[, c(which(colnames(somalia) %in% c("latitude", "longitude")), grep("Embedding", colnames(somalia)), which(colnames(somalia) == "duration"))]
    importance = sqrt(importance/c(2, 6, 1))
    data = apply(data, 2, function(x){x - mean(x)})
    data[, c(1, 2)] = data[, c(1, 2)]/sqrt(var(as.vector(as.matrix(data[, c(1, 2)])))) * importance[1]
    data[, 3:(dim(data)[2] - 1)] = data[, 3:(dim(data)[2] - 1)]/sqrt(var(as.vector(as.matrix(data[, 3:(dim(data)[2] - 1)])))) * importance[2]
    data[, dim(data)[2]] = data[, dim(data)[2]]/sqrt(var(data[, dim(data)[2]])) * importance[3]
    return(as.matrix(data))
}
mySpecc = function(dataNormalized, numClusters, filePath, forceRun = FALSE){
    if(!file.exists(filePath) | forceRun){
        ClustersObj = specc(dataNormalized, centers = numClusters)
        saveRDS(ClustersObj, filePath)
    }else{
        ClustersObj = readRDS(filePath)
    }
    return(ClustersObj)
}
calculateVariableContribution = function(dataNormalized, Clusters){
    dataGroups = split(as.data.frame(dataNormalized), factor(Clusters))
    withinMSE = rowSums(sapply(dataGroups, function(group){apply(group, 2, function(x){sum((x - mean(x))^2)})}))/(dim(somalia)[1] - numClusters)
    withinMSE = c(sum(withinMSE[1:2]), sum(withinMSE[-c(1, 2, length(withinMSE))]), withinMSE[length(withinMSE)])
    betweenMSE = rowSums((sapply(dataGroups, colMeans) - colMeans(dataNormalized))^2)/(numClusters - 1)
    betweenMSE = c(sum(betweenMSE[1:2]), sum(betweenMSE[-c(1, 2, length(betweenMSE))]), betweenMSE[length(betweenMSE)])
    importance = exp(-withinMSE/betweenMSE)
    contribution = importance/sum(importance)
    print(paste("Contributions of variables are: ", paste(round(contribution * 100, 2), collapse = "%, "), "%.", sep = ""))
    return(contribution)
}
getGroupsNear = function(data, Clusters, location = c(45.318161, 2.046934), returnGroupNum = 3, returnGroupSizeLargerThan = 1){
    # c(longitude = 45.318161, latitude = 2.046934) is the location of Mogadishu, data = somalia; Clusters = somaliaClusters
    distCenters2Location = sapply(split(as.data.frame(cbind(data$longitude, data$latitude)), Clusters), function(x){sum((colMeans(x) - location)^2)})
    distCenters2Location[table(Clusters) < returnGroupSizeLargerThan] = max(distCenters2Location) 
    return(as.numeric(names(sort(distCenters2Location)[1:returnGroupNum])))
}
getDataFromGroups = function(data, Clusters, groups){
    result = vector("list", length = length(groups))
    for(i in 1:length(groups)){
        gp = groups[i]
        result[[i]] = data.frame(actor1 = data$actor1[Clusters == gp], 
                                 actor2 = data$actor2[Clusters == gp],
                                 duration = data$duration[Clusters == gp])
    }
    return(result)
}
plotGroupsNear = function(data, Clusters, groups, color_vec = c("blue", "green", "red", "yellow", "purple", "orange", "black"), ...){
    plot(data$longitude, data$latitude, cex = 1, pch = 20, col = "grey", ...)
    for(i in 1:length(groups)){
        gp = groups[i]
        points(data$longitude[Clusters == gp], data$latitude[Clusters == gp], cex = 2, pch = 20, col = color_vec[i])
    }
}
